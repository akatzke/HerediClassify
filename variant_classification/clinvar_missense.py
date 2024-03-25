#!/usr/bin/env python3

from math import ceil
import logging
import pathlib
import pandas as pd

from Bio.Seq import Seq
from Bio.Data import IUPACData
from cyvcf2 import VCF
import pyensembl

from variant import VariantInfo, TranscriptInfo
from var_type import VARTYPE_GROUPS
from genotoscope_exon_skipping import (
    is_transcript_in_positive_strand,
)
from clinvar_utils import (
    ClinVar,
    ClinVar_Type,
    filter_gene,
    create_ClinVar,
    convert_vcf_gen_to_df,
    get_affected_transcript,
)
from custom_exceptions import No_transcript_with_var_type_found

logger = logging.getLogger("GenOtoScope_Classify.clinvar.missense")


def check_clinvar_missense(
    variant: VariantInfo, transcripts: list[TranscriptInfo], path_clinvar: pathlib.Path
) -> tuple[ClinVar, ClinVar]:
    """
    Check ClinVar for entries supporting pathogenicity of missense variant
    """
    # In case the variant is not an SNV, return empty ClinVar result
    if len(variant.var_obs) != 1 or len(variant.var_ref) != 1:
        logger.warning(
            "Variant is not a SNV. PS1/PM5 currently not implemented for delins."
        )
        ClinVar_same_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.SAME_AA_CHANGE)
        ClinVar_diff_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.DIFF_AA_CHANGE)
        return (ClinVar_same_aa, ClinVar_diff_aa)
    try:
        affected_transcript, ref_transcript = get_affected_transcript(
            transcripts, VARTYPE_GROUPS.MISSENSE
        )
    except No_transcript_with_var_type_found:
        logger.warning("No transcript with variant type missense found.")
        ClinVar_same_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.SAME_AA_CHANGE)
        ClinVar_diff_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.DIFF_AA_CHANGE)
        return (ClinVar_same_aa, ClinVar_diff_aa)

    var_codon_info = extract_var_codon_info(
        variant, ref_transcript, affected_transcript
    )
    clinvar_same_codon = extract_clinvar_entries_missense(
        path_clinvar,
        variant.chr,
        var_codon_info["genomic_pos"],
        var_codon_info["intersects_intron_at"],
    )
    clinvar_same_codon_aa: pd.DataFrame = clinvar_same_codon.apply(
        construct_clinvar_prot_change,
        var_codon_info=var_codon_info,
        axis=1,
    )
    # Filter out all ClinVar entries with termination codon
    if not clinvar_same_codon_aa.empty:
        clinvar_same_codon_no_ter = clinvar_same_codon_aa[
            clinvar_same_codon_aa.prot_alt != "Ter"
        ]
    else:
        clinvar_same_codon_no_ter = clinvar_same_codon_aa
    if not clinvar_same_codon_no_ter.empty:
        clinvar_same_codon_aa_filtered = filter_gene(
            clinvar_same_codon_no_ter, variant.gene_name
        )
        clinvar_same_aa_df = clinvar_same_codon_aa_filtered[
            clinvar_same_codon_aa_filtered.prot_alt == var_codon_info["prot_alt"]
        ]
        clinvar_diff_aa = clinvar_same_codon_aa_filtered[
            clinvar_same_codon_aa_filtered.prot_alt != var_codon_info["prot_alt"]
        ]
    else:
        clinvar_same_aa_df = pd.DataFrame()
        clinvar_diff_aa = pd.DataFrame()
    ClinVar_same_aa = create_ClinVar(clinvar_same_aa_df, ClinVar_Type.SAME_AA_CHANGE)

    ClinVar_diff_aa = create_ClinVar(clinvar_diff_aa, ClinVar_Type.DIFF_AA_CHANGE)

    return (ClinVar_same_aa, ClinVar_diff_aa)


def extract_var_codon_info(
    variant_info: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    transcript: TranscriptInfo,
) -> dict:
    """
    Construct variant codon information:
    genomic_pos: list[int]
    intersect_intron_at: int
    var_strand: str
    seq_ref: str
    prot_start: int
    prot_ref: str
    prot_alt str
    """
    logger.debug("Extract codon information from variant-affected genomic position")
    logger.debug(f"=== New transcript id = {transcript.transcript_id} ===")
    var_start = int(transcript.var_hgvs.pos.start.base)
    var_pos_in_codon = get_variant_position_in_codon(var_start)

    if is_transcript_in_positive_strand(ref_transcript):
        var_strand = "+"
    else:
        var_strand = "-"

    logger.debug(f"var_pos_in_codon = {var_pos_in_codon}")

    codon_genomic_pos = construct_codon_position(
        variant_info.genomic_start, var_pos_in_codon, var_strand
    )
    codon_coding_pos = construct_codon_position(var_start, var_pos_in_codon, "+")
    logger.debug(f"var_codon_genomic_pos: {codon_genomic_pos}")

    # correct variant codon genomic positions for 2 exon spanning codons
    (
        var_codon_genomic_pos_corrected,
        codon_intersects_intron_at,
    ) = normalize_codon_exonic_pos(
        ref_transcript, variant_info, codon_genomic_pos, var_strand
    )
    ### ### ### ### ### ### ### ### ###
    # create reference codon sequence #
    ### ### ### ### ### ### ### ### ###
    try:
        # use the coding sequence attribute
        codon_seq_ref = ref_transcript.coding_sequence[
            codon_coding_pos[0] - 1 : codon_coding_pos[-1]
        ]
    except ValueError:
        # solve value error by using the sequence attribute instead
        # Pyensembl returned ValueError for current transcript with id (seen for mitochondrial genes)
        codon_seq_ref = ref_transcript.sequence[
            codon_coding_pos[0] - 1 : codon_coding_pos[-1]
        ]
        logger.debug(
            f"Use sequence attribute instead of coding_sequence for transcript id: {transcript.transcript_id}"
        )

    codon_seq_obs = construct_obs_codon_seq(
        codon_seq_ref, var_pos_in_codon, variant_info.var_obs, var_strand
    )

    # create amino of protein for affected codon
    prot_var_start = ceil(var_start / 3)
    prot_ref = convert_codon_to_aa(codon_seq_ref)
    prot_alt = convert_codon_to_aa(codon_seq_obs)

    var_codon_info = {
        "genomic_pos": var_codon_genomic_pos_corrected,
        "intersects_intron_at": codon_intersects_intron_at,
        "strand": var_strand,
        "seq_ref": codon_seq_ref,
        "prot_start": prot_var_start,
        "prot_ref": prot_ref,
        "prot_alt": prot_alt,
    }
    logger.debug(f"Var codon info per transcript: {var_codon_info}")
    return var_codon_info


def get_variant_position_in_codon(var_start: int) -> int:
    """
    From cDNA position extrapolate position of variant in codon
    1-based
    """
    codon_index = var_start // 3
    var_pos_in_codon = var_start % 3

    if not var_start == codon_index * 3 + var_pos_in_codon:
        raise AssertionError(
            f"AssertionError: codon_index * 3 + var_pos_in_codon should be equal to codon index\n=> variant position: {var_start}"
        )
    if var_pos_in_codon == 0:
        var_pos_in_codon = 3
    logger.debug(
        f"Var start pos: {var_start} is at codon index: {codon_index} and variant is the codon position of {var_pos_in_codon}"
    )
    if not var_start >= var_pos_in_codon:
        raise AssertionError(
            f"Variant can't be at position {var_pos_in_codon} of codon and not be at least at the coding position {var_pos_in_codon}\n=> variant position: {var_start}"
        )
    return var_pos_in_codon


def construct_clinvar_prot_change(
    clinvar_rec: pd.Series, var_codon_info: dict
) -> pd.Series:
    """
    Construct protein change for clinvar record in dataframe

    For termination codon the return protein change will contain the 3-letter protein letter 'Ter'
    as discussed: https://varnomen.hgvs.org/recommendations/protein/variant/frameshift/
    """
    logger.debug(f"Construct protein change for clinvar record: {clinvar_rec}")
    # get position of clinvar nucleotide
    clinvar_rec["var_pos_idx"] = get_variant_position_in_codon_from_genomic_position(
        var_codon_info["strand"], var_codon_info["genomic_pos"], clinvar_rec.pos
    )
    logger.debug(
        f"ClinVar record is found in codon position: {clinvar_rec.var_pos_idx}"
    )
    ### ### ###
    # construct the reference coding sequence for clinvar record
    ### ### ###
    clinvar_ref_seq = list(var_codon_info["seq_ref"])
    logger.debug(f"clinvar ref: {clinvar_ref_seq}")
    ### ### ###
    # and then construct the alternate sequence
    ### ### ###
    clinvar_alt_seq = []
    corrected_obs_base = correct_observed_base_for_strand(
        var_codon_info["strand"], clinvar_rec.alt
    )
    for idx, nucl in enumerate(var_codon_info["seq_ref"]):
        if idx == clinvar_rec.var_pos_idx:
            clinvar_alt_seq.append(corrected_obs_base)
        else:
            clinvar_alt_seq.append(nucl)
    logger.debug(f"clinvar alt: {clinvar_alt_seq}")

    ### ### ###
    # translate these ref and alt sequence to protein edit
    ### ### ###
    clinvar_ref_translated = str(
        Seq(normalize_codon_nucleotides(clinvar_ref_seq)).translate()
    )
    clinvar_rec["prot_ref"] = "".join(convert_1to3_aa(clinvar_ref_translated))
    clinvar_alt_translated = str(
        Seq(normalize_codon_nucleotides(clinvar_alt_seq)).translate()
    )
    clinvar_rec["prot_alt"] = "".join(convert_1to3_aa(clinvar_alt_translated))
    logger.debug(
        f"ref_translated: {clinvar_ref_translated}, alt_translated: {clinvar_alt_translated}"
    )
    clinvar_rec["prot_start"] = var_codon_info["prot_start"]
    return clinvar_rec


def normalize_codon_nucleotides(nucl_per_codon_pos) -> Seq:
    """
    Normalize nucleotides registered for the three positions of a codon

    Parameters
    ----------
    nucl_per_codon_pos: list of str
        nucleotides per codon position
    """

    logger.debug(
        "Normalize nucleotides: {}, registered on one codon".format(nucl_per_codon_pos)
    )
    seq2translate = ""
    nucl_pos = (nucl for nucl in nucl_per_codon_pos)
    while len(seq2translate) < 3:
        seq2translate += next(nucl_pos)

    logger.debug("Normalized pos: {}".format(seq2translate))
    logger.debug("Padded: {}".format(pad_seq(seq2translate)))
    return pad_seq(seq2translate)


def convert_1to3_aa(amino_acid: str) -> str:
    """
    Convert 1 letter amino acid genotoscope to 3 letter equivalent
    For termination codon the return protein change will contain the 3-letter protein letter 'Ter'
    as discussed: https://varnomen.hgvs.org/recommendations/protein/variant/frameshift/
    """
    try:
        amino_acid_3let = IUPACData.protein_letters_1to3[amino_acid]
    except KeyError:
        amino_acid_3let = "Ter"
    return amino_acid_3let


def extract_clinvar_entries_missense(
    path_clinvar: pathlib.Path,
    chrom: str,
    genomic_positions: list[int],
    codon_intersects_intron: bool,
) -> pd.DataFrame:
    """
    Extract matching ClinVar entries for missense variants
    """
    clinvar = VCF(path_clinvar)
    if not codon_intersects_intron:
        clinvar_same_codon = clinvar(
            f"{chrom}:{genomic_positions[0]}-{genomic_positions[2]}"
        )
        clinvar_same_codon_df = convert_vcf_gen_to_df(clinvar_same_codon)
    else:
        clinvar_same_codon_df = pd.DataFrame()
        for position in genomic_positions:
            clinvar_one_nucleotide = clinvar(f"{chrom}:{position}-{position}")
            clinvar_df = convert_vcf_gen_to_df(clinvar_one_nucleotide)
            clinvar_same_codon_df = pd.concat([clinvar_same_codon_df, clinvar_df])
    return clinvar_same_codon_df


def construct_codon_position(
    var_start: int, var_pos_in_codon: int, strand: str
) -> list[int]:
    """
    Construct genomic position based on variant genomic start and variant codon position
    Construct cDNA position of variant
    """
    if var_pos_in_codon == 2:
        codon_pos = [
            var_start - 1,
            var_start,
            var_start + 1,
        ]
    elif var_pos_in_codon == 1:
        if strand == "-":
            codon_pos = [
                var_start - 2,
                var_start - 1,
                var_start,
            ]
        else:
            codon_pos = [
                var_start,
                var_start + 1,
                var_start + 2,
            ]
    elif var_pos_in_codon == 3:
        if strand == "-":
            codon_pos = [
                var_start,
                var_start + 1,
                var_start + 2,
            ]
        else:
            codon_pos = [
                var_start - 2,
                var_start - 1,
                var_start,
            ]
    else:
        assert (
            var_pos_in_codon <= 3 and var_pos_in_codon >= 1
        ), f"Assertion error: Variant codon postition should be between 1 and 4 \n Variant codon position is {var_pos_in_codon}"
    return codon_pos


def normalize_codon_exonic_pos(
    ref_transcript: pyensembl.transcript.Transcript,
    variant_info: VariantInfo,
    codon_genomic_positions: list[int],
    transcript_strand: str,
) -> tuple:
    """
    Correct codon position after investigating intersection with an intron
    Returns
    -------
    codon_pos_corrected : list of int
        corrected codon positions
    codon_intersects_intron : bool
    """
    logger.debug("Normalize codon in exonic positions")
    ### ### ###
    # create coding sequence ranges respecting transcript's strand direction
    ### ### ###
    if transcript_strand == "+":
        transcript_strand = +1
        codon_start, codon_end = codon_genomic_positions[0], codon_genomic_positions[2]
    else:
        transcript_strand = -1
        codon_start, codon_end = codon_genomic_positions[2], codon_genomic_positions[0]
    codon_middle = codon_genomic_positions[1]

    # for coding_range in transcript.coding_sequence_position_ranges:
    codon_intersects_intron_at = -1
    coding_seq_positions = []
    for exon in ref_transcript.exons:
        if transcript_strand == +1:
            coding_seq_positions.append(
                [exon.to_dict()["start"], exon.to_dict()["end"]]
            )
        else:
            coding_seq_positions.append(
                [exon.to_dict()["end"], exon.to_dict()["start"]]
            )

    ### ### ### ### ### ### ### ### ### ### ### ### ###
    # find at which exon (or between which two exons) #
    # the codon lies in                               #
    ### ### ### ### ### ### ### ### ### ### ### ### ###
    codon_pos_corrected_tmp = []
    for cod_seq_idx, coding_seq_range in enumerate(coding_seq_positions):
        coding_seq_start, coding_seq_end = coding_seq_range[0], coding_seq_range[1]
        normalized_exon_interval = range(
            coding_seq_start, coding_seq_end + transcript_strand, transcript_strand
        )

        # logger.debug(f"Coding start: {coding_seq_start}, end: {coding_seq_end}, range: {normalized_exon_interval}")
        # if codon_start >= coding_seq_start and codon_end <= coding_seq_end:
        if (
            codon_start in normalized_exon_interval
            and codon_end in normalized_exon_interval
        ):
            logger.debug("codon in exon")
            codon_pos_corrected = codon_genomic_positions
            codon_intersects_intron_at = -1
            break
        elif codon_start in normalized_exon_interval:
            logger.debug("start in exon")
            if codon_start == coding_seq_end:
                logger.debug("Codon start is the last base of exon")
                if transcript_strand == +1:
                    codon_pos_corrected_tmp.append([codon_start])
                    codon_pos_corrected_tmp.append(
                        [
                            coding_seq_positions[cod_seq_idx + 1][0],
                            coding_seq_positions[cod_seq_idx + 1][0] + 1,
                        ]
                    )
                else:
                    # codon start is the one before the last base of the exon
                    codon_pos_corrected_tmp.append(
                        [
                            coding_seq_positions[cod_seq_idx + 1][0] - 1,
                            coding_seq_positions[cod_seq_idx + 1][0],
                        ]
                    )
                    codon_pos_corrected_tmp.append([codon_start])
                # intersects on the 2nd position (0-index = 1)
                codon_intersects_intron_at = 1
                codon_pos_corrected = sum(codon_pos_corrected_tmp, [])
            elif codon_start == coding_seq_end - 1 * transcript_strand:
                logger.debug("Codon start is the penultimate base of exon")
                if transcript_strand == +1:
                    codon_pos_corrected_tmp.append([codon_start, codon_middle])
                    codon_pos_corrected_tmp.append(
                        [coding_seq_positions[cod_seq_idx + 1][0]]
                    )
                else:
                    codon_pos_corrected_tmp.append(
                        [coding_seq_positions[cod_seq_idx + 1][0]]
                    )
                    codon_pos_corrected_tmp.append([codon_middle, codon_start])
                # intersects on the 3rd position (0-index = 2)
                codon_intersects_intron_at = 2
                codon_pos_corrected = sum(codon_pos_corrected_tmp, [])
            break
        elif (
            codon_middle in normalized_exon_interval
            and codon_end in normalized_exon_interval
        ):
            logger.debug("middle,end in exon")
            if transcript_strand == +1:
                codon_pos_corrected_tmp.append(
                    [coding_seq_positions[cod_seq_idx - 1][1]]
                )
                codon_pos_corrected_tmp.append([codon_middle, codon_end])
            else:
                codon_pos_corrected_tmp.append([codon_end, codon_middle])
                codon_pos_corrected_tmp.append(
                    [coding_seq_positions[cod_seq_idx - 1][1]]
                )
            # genomic position of codon intersects 1st position (0-index = 0)
            codon_intersects_intron_at = 0
            codon_pos_corrected = sum(codon_pos_corrected_tmp, [])
            break
        elif codon_end in normalized_exon_interval:
            logger.debug("end in exon")
            if transcript_strand == +1:
                codon_pos_corrected_tmp.append(
                    [
                        coding_seq_positions[cod_seq_idx - 1][1] - 1,
                        coding_seq_positions[cod_seq_idx - 1][1],
                    ]
                )
                codon_pos_corrected_tmp.append([codon_end])
            else:
                codon_pos_corrected_tmp.append([codon_end])
                codon_pos_corrected_tmp.append(
                    [
                        coding_seq_positions[cod_seq_idx - 1][1],
                        coding_seq_positions[cod_seq_idx - 1][1] + 1,
                    ]
                )
            # genomic position of codon intersects 1st position (0-index = 0)
            codon_intersects_intron_at = 0
            codon_pos_corrected = sum(codon_pos_corrected_tmp, [])
            break
        else:
            continue
    logger.debug(f"normalized positions: {codon_pos_corrected}")
    logger.debug(f"codon intersects intron at: {codon_intersects_intron_at}")
    try:
        assert len(codon_pos_corrected) == 3
    except AssertionError:
        logger.error(
            f"AssertionError: Codon does not intersects intron, thus corrected codon positions should be 3 \n=> variant position: {variant_info.to_string()}",
            exc_info=True,
        )
    if codon_intersects_intron_at == -1:
        # assert corrected positions are sorted
        try:
            assert sorted(codon_pos_corrected) == codon_pos_corrected
        except AssertionError:
            logger.error(
                f"AssertionError: Corrected position of codon are not sorted: {codon_pos_corrected}\n=> variant position: {variant_info.to_string()}",
                exc_info=True,
            )
    else:
        # check that codon positions are increasing
        start_pos = 0
        for codon_position in codon_pos_corrected:
            try:
                assert codon_position > start_pos
            except AssertionError:
                logger.error(
                    f"AssertionError: Codon positions are not increasing: {codon_pos_corrected}\n=> variant position: {variant_info.to_string}",
                    exc_info=True,
                )
                start_pos = codon_position
    if codon_intersects_intron_at == -1:
        codon_intersects_intron = False
    else:
        codon_intersects_intron = True
    return codon_pos_corrected, codon_intersects_intron


def get_variant_position_in_codon_from_genomic_position(
    var_strand: str, codon_pos_genomic: list[int], var_genomic_pos: int
) -> int:
    """
    From variant genomic position and genomic positions of codon get var_pos_in_codon
    0-based index
    """
    if var_strand == "+":
        var_pos_in_codon = codon_pos_genomic.index(int(var_genomic_pos))
        return var_pos_in_codon
    else:
        convert_index = {0: 2, 1: 1, 2: 0}
        var_pos_in_codon_pos_strand = codon_pos_genomic.index(int(var_genomic_pos))
        var_pos_in_codon_neg_strand = convert_index[var_pos_in_codon_pos_strand]
        return var_pos_in_codon_neg_strand


def pad_seq(sequence: Seq) -> Seq:
    """
    Pad sequence to multiple of 3 with N
    Credits: https://stackoverflow.com/questions/53894575/how-can-i-fix-this-error-biopythonwarning-partial-codon-lensequence-not-a
    """

    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + "N" * (3 - remainder)


def construct_obs_codon_seq(
    codon_seq_ref: str, var_pos_in_codon: int, obs_base: str, strand: str
) -> str:
    """
    Construct observed codon sequence
    """
    corrected_obs_base = correct_observed_base_for_strand(strand, obs_base)
    if var_pos_in_codon == 3:
        codon_seq_obs = [char for char in codon_seq_ref[0:2]] + [corrected_obs_base]
    elif var_pos_in_codon == 1:
        codon_seq_obs = [corrected_obs_base] + [char for char in codon_seq_ref[1:3]]
    elif var_pos_in_codon == 2:
        codon_seq_obs = [codon_seq_ref[0], corrected_obs_base, codon_seq_ref[2]]
    codon_seq_obs = "".join(codon_seq_obs)

    # assert that the codon size is multiple of 3
    try:
        len(codon_seq_ref) % 3 == 0 and len(codon_seq_obs) % 3 == 0
    except AssertionError:
        logger.error(
            f"The codon sequence for reference or observed sequence is not multiple of 3",
            exc_info=True,
        )
    return codon_seq_obs


def convert_codon_to_aa(codon: str) -> str:
    """
    Get 3 letter amino acid code corresponding to given codon
    """
    aa_1let = str(Seq(codon).translate())
    aa_3let = convert_1to3_aa(aa_1let)
    return aa_3let


def correct_observed_base_for_strand(strand: str, base: str) -> str:
    """
    Depending on which strand the transcript is located on, correct observed base
    """
    comp_bases = {"C": "G", "G": "C", "A": "T", "T": "A"}
    if strand == "+":
        corrected_base = base
    else:
        corrected_base = ""
        for entry in base:
            corrected_base += comp_bases[entry]
    return corrected_base
