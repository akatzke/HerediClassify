#!/usr/bin/env python3

import logging
import pathlib
import pandas as pd
from typing import Generator
from math import ceil
from cyvcf2 import VCF
from Bio.Data import IUPACData
from Bio.Seq import Seq
import pyensembl
import hgvs.parser
import hgvs.posedit
from refactoring.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
)

from refactoring.variant import VariantInfo, TranscriptInfo
from refactoring.genotoscope_exon_skipping import is_transcript_in_positive_strand

hgvs_parser = hgvs.parser.Parser()

logger = logging.getLogger("GenOtoScope_Classify.PVS1.assess_NMD")


def check_clinvar_splicing(
    variant: VariantInfo,
    transcripts: list[TranscriptInfo_intronic],
    path_clinvar: pathlib.Path,
) -> tuple:
    """
    Check ClinVar for entries supporting pathogenicity of splice site
    """
    clinvar = VCF(path_clinvar)
    clinvar_same_pos = clinvar(
        f"{variant.chr}:{variant.genomic_start}-{variant.genomic_end}"
    )
    clinvar_same_pos_df = convert_vcf_gen_to_df(clinvar_same_pos)
    if not clinvar_same_pos_df.empty:
        max_classification, affected_ID = get_highest_classification(clinvar_same_pos)
        return True, True
        # return ClinVar_splice(same_nucleotide_change_pathogenic=True, matching_clinvar_entries=clinvar_same_pos_filtered, clinvar_pos_fil)
    else:
        (start_splice_site, end_splice_site) = find_corresponding_splice_site(variant)
        clinvar_splice_site = clinvar(
            f"{variant.chr}:{start_splice_site}-{end_splice_site}"
        )
        clinvar_splice_site_df = convert_vcf_gen_to_df(clinvar_splice_site)
        if clinvar_splice_site:
            return True, False
        else:
            return False, False


def check_clinvar_missense(
    variant: VariantInfo, transcripts: list[TranscriptInfo_exonic], path_clinvar: str
):
    """
    Check ClinVar for entries supporting pathogenicity of missense variant
    """
    var_codon_info = extract_var_codon_info(variant, transcripts)


def convert_vcf_gen_to_df(vcf_generator: Generator) -> pd.DataFrame:
    """
    Covnerts cyvcf generator into a pd.DataFrame
    """
    names = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]
    df = pd.DataFrame(columns=names)
    for entry in vcf_generator:
        clinvar_split_str = str(entry).split("\t")
        clinvar_dict = dict(zip(names, clinvar_split_str))
        df = pd.concat([df, pd.DataFrame([clinvar_dict])], axis=0, ignore_index=True)
    df_format_info = format_info(df)
    return df_format_info


def format_info(data: pd.DataFrame) -> pd.DataFrame:
    """
    Format the info column
    """
    info = data["info"]
    info_split = [entry.split("=") for entry in info]
    info_split = [entry.split("\n")[0].split(";") for entry in info]
    processed_info = [[item.split("=") for item in entry] for entry in info_split]
    dict_info = [
        {item[0]: item[1] for item in entry if len(item) > 1}
        for entry in processed_info
    ]
    return pd.concat([data, pd.DataFrame(dict_info)], axis=1)


def check_clinvar_region(chr: int, start: int, end: int, path_clinvar: str):
    clinvar_df = VCF(path_clinvar)
    clinvar_region = clinvar_df(f"{chr}:{start}-{end}")
    return clinvar_region


def find_corresponding_splice_site(variant: VariantInfo) -> tuple:
    return (variant.genomic_start, variant.genomic_end)


def get_highest_classification(clinvars) -> str:
    return clinvars.ID


def extract_var_codon_info(
    variant_info: VariantInfo, transcripts: list[TranscriptInfo_exonic]
) -> dict[str, dict[str, str]]:
    """
    Get location of variant affected codon in each transcript
    """
    logger.debug("Extract codon information from variant-affected genomic position")
    transcripts_var_codon_info = {}
    for transcript in transcripts:
        logger.debug(f"=== New transcript id = {transcript.transcript_id} ===")
        ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            transcript.transcript_id
        )
        # find at which codon is the variant located
        # based on the variant location in the coding cDNA
        var_start = int(transcript.var_hgvs.pos.start.base)
        codon_index = var_start // 3
        var_codon_pos = var_start % 3

        try:
            var_start == codon_index * 3 + var_codon_pos
        except AssertionError:
            logger.error(
                f"AssertionError: codon_index * 3 + var_codon_pos should be equal to codon index\n=> variant position: {variant_info.to_string()}",
                exc_info=True,
            )
        if var_codon_pos == 0:
            var_codon_pos = 3
        logger.debug(
            f"Var start pos: {var_start} is at codon index: {codon_index} and variant is the codon position of {var_codon_pos}"
        )
        var_codon_genomic_pos, var_codon_coding_pos = [], []
        logger.debug(f"var_codon_pos = {var_codon_pos}")

        # understand variant edit type
        # create genomic position and amino acid per codon
        if var_codon_pos == 3:  # variant is at the third (last) position of a codon
            try:
                var_start >= 3
            except AssertionError:
                logger.error(
                    f"Variant can't be at the last position of codon and not be at least at the coding position 3\n=> variant position: {variant_info.to_string()}",
                    exc_info=True,
                )
            var_codon_genomic_pos = [
                variant_info.genomic_start - 2,
                variant_info.genomic_start - 1,
                variant_info.genomic_start,
            ]
            var_codon_coding_pos = [var_start - 2, var_start - 1, var_start]
        elif var_codon_pos == 1:  # variant is at the first position of a codon
            var_codon_genomic_pos = [
                variant_info.genomic_start,
                variant_info.genomic_start + 1,
                variant_info.genomic_start + 2,
            ]
            var_codon_coding_pos = [var_start, var_start + 1, var_start + 2]
        elif var_codon_pos == 2:  # variant is at the second position of a codon
            try:
                var_start >= 2
            except AssertionError:
                logger.error(
                    f"Variant can't be at the second position of a codon\n=> variant position: {variant_info.to_string()}",
                    exc_info=True,
                )

            var_codon_genomic_pos = [
                variant_info.genomic_start - 1,
                variant_info.genomic_start,
                variant_info.genomic_start + 1,
            ]
            var_codon_coding_pos = [var_start - 1, var_start, var_start + 1]
        logger.debug(f"var_codon_genomic_pos: {var_codon_genomic_pos}")

        # get the transcript strand
        if is_transcript_in_positive_strand(ref_transcript):
            var_strand = "+"
        else:
            var_strand = "-"

        # correct variant codon genomic positions
        (
            var_codon_genomic_pos_corrected,
            codon_intersects_intron_at,
        ) = normalize_codon_exonic_pos(
            ref_transcript, variant_info, var_codon_genomic_pos, var_strand
        )

        ### ### ### ### ### ### ### ### ###
        # create reference codon sequence #
        ### ### ### ### ### ### ### ### ###
        try:
            # use the coding sequence attribute
            codon_seq_ref = ref_transcript.coding_sequence[
                var_codon_coding_pos[0] - 1 : var_codon_coding_pos[-1]
            ]
        except ValueError:
            # solve value error by using the sequence attribute instead
            # Pyensembl returned ValueError for current transcript with id (seen for mitochondrial genes)
            codon_seq_ref = ref_transcript.sequence[
                var_codon_coding_pos[0] - 1 : var_codon_coding_pos[-1]
            ]
            logger.debug(
                f"Use sequence attribute instead of coding_sequence for transcript id: {transcript_info['transcript_id']}"
            )

        ### ### ### ### ### ### ### ### ###
        # create observed codon sequence  #
        ### ### ### ### ### ### ### ### ###
        if var_codon_pos == 3:  # variant on the last position of the codon
            codon_seq_obs = [char for char in codon_seq_ref[0:2]] + [
                variant_info.var_obs
            ]
        elif var_codon_pos == 1:  # variant on the first position of the codon
            codon_seq_obs = [variant_info.var_obs] + [
                char for char in codon_seq_ref[1:3]
            ]
        elif var_codon_pos == 2:  # variant on the middle (second) position of the codon
            codon_seq_obs = [codon_seq_ref[0], variant_info.var_obs, codon_seq_ref[2]]
            codon_seq_obs = "".join(codon_seq_obs)

        # assert that the codon size is multiple of 3
        try:
            len(codon_seq_ref) % 3 == 0 and len(codon_seq_obs) % 3 == 0
        except AssertionError:
            logger.error(
                f"The codon sequence for reference or observed sequence is not multiple of 3\n=> variant position: {variant_info.to_string()}",
                exc_info=True,
            )

        # create amino of protein for affected codon
        prot_var_start = ceil(var_start / 3)
        codon_amino_ref, codon_amino_obs = [], []
        codon_amino_ref.append(str(Seq(codon_seq_ref).translate()))
        codon_amino_obs.append(str(Seq(codon_seq_obs).translate()))
        codon_amino_ref.append("".join(convert_1to3_aa(codon_amino_ref[0])))
        codon_amino_obs.append("".join(convert_1to3_aa(codon_amino_obs[0])))

        transcripts_var_codon_info[transcript.transcript_id] = {
            "var_start": var_start,
            "genomic_pos": var_codon_genomic_pos_corrected,
            "coding_pos": var_codon_coding_pos,
            "intersects_intron_at": codon_intersects_intron_at,
            "strand": var_strand,
            "seq_ref": codon_seq_ref,
            "seq_obs": codon_seq_obs,
            "prot_start": prot_var_start,
            "amino_ref": codon_amino_ref,
            "amino_obs": codon_amino_obs,
        }
        logger.debug(f"Var codon info per transcript: {transcripts_var_codon_info}")
    return transcripts_var_codon_info


def normalize_codon_exonic_pos(
    transcript, variant_info, codon_genomic_positions, transcript_strand
):
    """
    Correct codon position after investigating intersection with an intron

    Parameters
    ----------
    transcript : pyensembl.transcript
        transcript ensembl object
    variant_info : VariantInfo
        basic variant information
    codon_genomic_positions : list of int
        initial codon genomic positions
    transcript_strand : str
        strand of transcript

    Returns
    -------
    codon_pos_corrected : list of list of int
        corrected codon positions
    codon_intersects_intron_at : int
        position that codon intersects with an intron (0-index)
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
    coding_seq_positions = []
    for exon in transcript.exons:
        if transcript_strand == +1:
            coding_seq_positions.append(
                [exon.to_dict()["start"], exon.to_dict()["end"]]
            )
        else:
            coding_seq_positions.append(
                [exon.to_dict()["end"], exon.to_dict()["start"]]
            )
    num_coding_sequences = len(coding_seq_positions)

    ### ### ### ### ### ### ### ### ### ### ### ### ###
    # find at which exon (or between which two exons) #
    # the codon lies in                               #
    ### ### ### ### ### ### ### ### ### ### ### ### ###
    codon_pos_corrected = []
    codon_intersects_intron_at = -1
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
            # codon inside coding seq =>
            # ( --- coding_seq --- )
            #      [ -codon- ]
            codon_pos_corrected = codon_genomic_positions
            codon_intersects_intron_at = -1
            break
        elif codon_start in normalized_exon_interval:
            logger.debug(
                "start in exon"
            )  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            # codon may intersect with current coding seq range in the coding seq range's right side
            # ( --- current_coding_seq --- ) ( -intron- ) ...
            #                             [ -codon- ]
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            if codon_start == coding_seq_end:
                logger.debug("Codon start is the last base of exon")
                if transcript_strand == +1:
                    codon_pos_corrected.append([codon_start])
                    codon_pos_corrected.append(
                        [
                            coding_seq_positions[cod_seq_idx + 1][0],
                            coding_seq_positions[cod_seq_idx + 1][0] + 1,
                        ]
                    )
                else:
                    # codon start is the one before the last base of the exon
                    codon_pos_corrected.append(
                        [
                            coding_seq_positions[cod_seq_idx + 1][0] - 1,
                            coding_seq_positions[cod_seq_idx + 1][0],
                        ]
                    )
                    codon_pos_corrected.append([codon_start])
                # intersects on the 2nd position (0-index = 1)
                codon_intersects_intron_at = 1
            elif codon_start == coding_seq_end - 1 * transcript_strand:
                logger.debug("Codon start is the penultimate base of exon")
                if transcript_strand == +1:
                    codon_pos_corrected.append([codon_start, codon_middle])
                    codon_pos_corrected.append(
                        [coding_seq_positions[cod_seq_idx + 1][0]]
                    )
                else:
                    codon_pos_corrected.append(
                        [coding_seq_positions[cod_seq_idx + 1][0]]
                    )
                    codon_pos_corrected.append([codon_middle, codon_start])
                # intersects on the 3rd position (0-index = 2)
                codon_intersects_intron_at = 2
            break
        elif (
            codon_middle in normalized_exon_interval
            and codon_end in normalized_exon_interval
        ):
            logger.debug("middle,end in exon")
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            # 	codon middle and end position intersects current exon
            # (--- previous_coding_seq ---) ( -intron- )   (--- current_coding_seq ---)
            #                       start               middle,end
            #                       [-------------------codon-----]
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            if transcript_strand == +1:
                codon_pos_corrected.append([coding_seq_positions[cod_seq_idx - 1][1]])
                codon_pos_corrected.append([codon_middle, codon_end])
            else:
                codon_pos_corrected.append([codon_end, codon_middle])
                codon_pos_corrected.append([coding_seq_positions[cod_seq_idx - 1][1]])
            # genomic position of codon intersects 1st position (0-index = 0)
            codon_intersects_intron_at = 0
            break
        elif codon_end in normalized_exon_interval:
            logger.debug("end in exon")
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            # codon end position intersects exon
            # (--- previous_coding_seq ---) ( -intron- )   (--- current_coding_seq ---) #                  start,middle                 end
            #                    [--------------------codon----]
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            if transcript_strand == +1:
                codon_pos_corrected.append(
                    [
                        coding_seq_positions[cod_seq_idx - 1][1] - 1,
                        coding_seq_positions[cod_seq_idx - 1][1],
                    ]
                )
                codon_pos_corrected.append([codon_end])
            else:
                codon_pos_corrected.append([codon_end])
                codon_pos_corrected.append(
                    [
                        coding_seq_positions[cod_seq_idx - 1][1],
                        coding_seq_positions[cod_seq_idx - 1][1] + 1,
                    ]
                )
            # genomic position of codon intersects 1st position (0-index = 0)
            codon_intersects_intron_at = 0
            break
    logger.debug(f"normalized positions: {codon_pos_corrected}")
    logger.debug(f"codon intersects intron at: {codon_intersects_intron_at}")
    if codon_intersects_intron_at == -1:
        try:
            assert len(codon_pos_corrected) == 3
        except AssertionError:
            logger.error(
                f"AssertionError: Codon does not intersects intron, thus corrected codon positions should be 3 \n=> variant position: {variant_info.to_string()}",
                exc_info=True,
            )
        # assert corrected positions are sorted
        try:
            assert sorted(codon_pos_corrected) == codon_pos_corrected
        except AssertionError:
            logger.error(
                f"AssertionError: Corrected position of codon are not sorted: {codon_pos_corrected}\n=> variant position: {variant_info.to_string()}",
                exc_info=True,
            )
    else:
        try:
            assert len(codon_pos_corrected) == 2
        except AssertionError:
            logger.error(
                f"AssertionError: Codon intersect intron at {codon_intersects_intron_at}, thus codon positions should two lists \n=> variant position: {variant_info.to_string()}",
                exc_info=True,
            )

        # check that codon positions are increasing
        start_pos = 0
        for codon_positions in codon_pos_corrected:
            for pos in codon_positions:
                try:
                    assert pos > start_pos
                except AssertionError:
                    logger.error(
                        f"AssertionError: Codon positions are not increasing: {codon_pos_corrected}\n=> variant position: {variant_info.to_string}",
                        exc_info=True,
                    )
                    start_pos = pos
    return codon_pos_corrected, codon_intersects_intron_at


def convert_1to3_aa(amino_acids: list[str]) -> list[str]:
    """
    Convert 1 letter amino acid genotoscope to 3 letter equivalent
    For termination codon the return protein change will contain the 3-letter protein letter 'Ter'
    as discussed: https://varnomen.hgvs.org/recommendations/protein/variant/frameshift/
    """

    aa_3code = []
    for aa in amino_acids:
        try:
            letter3 = IUPACData.protein_letters_1to3[aa]
        except KeyError:
            letter3 = "Ter"
        aa_3code.append(letter3)
    return aa_3code


def extract_clinvar_records(transcripts_var_codon_info, variant_info):
    """
    Extract ClinVar records matching codon genomic positions per transcript
    Each extracted clinvar record is a dictionary of
    ID -> clinvar id
    CLNDISDB -> clinvar disease db
    CLNREVSTAT -> converted gold stars from clinvar review status (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
    CLNSIG -> clinvar clinical significance

    Parameters
    ----------
    transcripts_var_codon_info : dict of str : dict of str : str
        two level dictionary, map each transcript id to codon information
    variant_info : VariantInfo
        current variant info object

    Returns
    -------
    dict of str : list of dict of str: str
        dictionary mapping transcript id to all clinvar record found for the specific codon of this transcript
    """

    logger.debug(
        "Extract ClinVar records that are found at the same codon as input variant"
    )
    codon_clinvars = {}
    for transcript_id, transcript_codon_info in transcripts_var_codon_info.items():
        logger.debug(f" === Transcript {transcript_id} === ")
        chr = variant_info.chrom.split("chr")[1]
        transcript = ensembl_data.transcript_by_id(transcript_id)
        transcript_strand = transcript.exons[0].to_dict()["strand"]
        ### ### ### ###
        # for each genomic position,
        # first search the cached clinvars
        # if no hit, then perform vcf fetch
        # vcf fetch uses 0-index and open upper limit (stop position)
        ### ### ### ###
        if transcript_codon_info["intersects_intron_at"] == -1:
            logger.debug("Codon does not intersect intron")
            # codon does not intersect intron
            start, end = (
                transcript_codon_info["genomic_pos"][0] - 1,
                transcript_codon_info["genomic_pos"][2],
            )
            matched_clinvars = []
            for pos in range(start, end):
                case_hit, cached_clinvars = match_pos2cache(chr, pos, transcript_strand)
                if not case_hit:
                    extracted_clinvars = add_info2extracted_clinvars(
                        list(vcf_reader.fetch(chr, pos, pos + 1))
                    )
                    logger.debug(f"Extracted clinvars: {extracted_clinvars}")
                    if len(extracted_clinvars) > 0:
                        # cache: save strand and quality filter ClinVars
                        # uniq_clinvars = uniq_id_filter(extracted_clinvars)
                        strand_filtered_clinvars = strand_filter(
                            extracted_clinvars, transcript_strand
                        )
                        quality_filtered_clinvars = quality_filter(
                            strand_filtered_clinvars, min_review_stars
                        )
                        unique_clinvars = uniq_new_filter(quality_filtered_clinvars)
                        cache_clinvars(chr, pos, transcript_strand, unique_clinvars)
                        matched_clinvars = matched_clinvars + unique_clinvars
                    else:
                        # cache: no ClinVar matching this chr,pos,strand
                        cache_clinvars(chr, pos, transcript_strand, [])
                        matched_clinvars = matched_clinvars + []
                else:
                    # use cached ClinVar entries
                    try:
                        assert len(cached_clinvars) >= 0
                    except AssertionError:
                        logger.error(
                            f"Number of extracted clinvars should be >= 0\n=> variant position: {variant_info.to_string()}",
                            exc_info=True,
                        )
                    matched_clinvars = matched_clinvars + cached_clinvars
        else:
            ### ### ###
            # aggregate clinvar entries for the 'separated' genomic positions of the codon
            ### ### ###
            logger.debug("Variant on two different exons")
            try:
                len(transcript_codon_info["genomic_pos"]) == 2
            except AssertionError:
                logger.error(
                    f"Corrected codon genomic position should be contained into two lists\n=> variant position: {variant_info.to_string}",
                    exc_info=True,
                )
            matched_clinvars = []
            for codon_genomic_range in transcript_codon_info["genomic_pos"]:
                if (
                    len(codon_genomic_range) == 2
                ):  # two positions are contained in the current codon partition
                    start, end = codon_genomic_range[0] - 1, codon_genomic_range[1]
                else:  # one position is contained in the current codon partition
                    start, end = codon_genomic_range[0] - 1, codon_genomic_range[0]
                for pos in range(start, end):
                    case_hit, cached_clinvars = match_pos2cache(
                        chr, pos, transcript_strand
                    )
                    if not case_hit:
                        extracted_clinvars = add_info2extracted_clinvars(
                            list(vcf_reader.fetch(chr, pos, pos + 1))
                        )
                        logger.debug(f"Extracted clinvars: {extracted_clinvars}")
                        if len(extracted_clinvars) > 0:
                            # cache: save strand and quality filter ClinVars
                            strand_filtered_clinvars = strand_filter(
                                extracted_clinvars, transcript_strand
                            )
                            quality_filtered_clinvars = quality_filter(
                                strand_filtered_clinvars, min_review_stars
                            )
                            unique_clinvars = uniq_new_filter(quality_filtered_clinvars)
                            cache_clinvars(chr, pos, transcript_strand, unique_clinvars)
                            matched_clinvars = matched_clinvars + unique_clinvars
                        else:
                            # cache: no ClinVar matching this chr,pos,strand
                            cache_clinvars(chr, pos, transcript_strand, [])
                            matched_clinvars = matched_clinvars + []
                    else:
                        # use cached ClinVars
                        try:
                            assert len(cached_clinvars) >= 0
                        except AssertionError:
                            logger.error(
                                f"Number of extracted clinvars should be >= 0\n=> variant position: {variant_info.to_string()}",
                                exc_info=True,
                            )
                        matched_clinvars = matched_clinvars + cached_clinvars
        if len(matched_clinvars) > 0:
            codon_clinvars[transcript_id] = uniq_new_filter(matched_clinvars)
    return codon_clinvars


def match_pos2cache(chr: str, pos: int, strand: int) -> tuple:
    """
    Match position to ClinVar cache

    Returns
    -------
    bool
        variant position found in cache (True), otherwise False
    list of dict of str: str
        cached ClinVars (or None)
    """

    logger.debug(f"Search pos: chr {chr}, {pos}, {strand} in cache")
    cache_hit, cached_clinvars = False, None
    if chr in clinvars_cache:
        if pos in clinvars_cache[chr]:
            if strand in clinvars_cache[chr][pos]:
                cache_hit = True
                cached_clinvars = clinvars_cache[chr][pos][strand]
    return cache_hit, cached_clinvars


def add_info2extracted_clinvars(clinvar_records):
    """
    Add information for extracted clinvar records
    Each extracted clinvar record is a dictionary of
    ID -> clinvar id
    CLNDISDB -> clinvar disease db
    CLNREVSTAT -> converted gold stars from clinvar review status (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
    CLNSIG -> clinvar clinical significance

    Parameters
    ----------
    clinvar_records : list of vcf.model._Record
        input clinvar records
    Returns
    -------
    list of dict
        clinvar records with needed information from vcf
    """

    logger.debug("Parse needed information for extracted ClinVar records")

    clinvars_info, uniq_ids = [], []
    for clinvar_rec in clinvar_records:
        if "CLNSIG" not in list(clinvar_rec.INFO.keys()):
            # if significance is not specified, pass
            continue
        # parse disease db
        if "CLNDISDB" in clinvar_rec.INFO:
            if clinvar_rec.INFO["CLNDISDB"][0]:
                rec_dis_db = ",".join(clinvar_rec.INFO["CLNDISDB"])
            else:
                rec_dis_db = "not_specified"
        else:
            rec_dis_db = "not_specified"
        clinvar_dict = {
            "id": clinvar_rec.ID,
            "pos": int(clinvar_rec.POS),
            "ref": str(clinvar_rec.REF),
            "alt": ",".join([str(alt) for alt in clinvar_rec.ALT]),
            "CLNDISDB": rec_dis_db,
            "CLNREVSTAT": convert_review_status2stars(
                clinvar_stars_df, star_status2int, clinvar_rec.INFO["CLNREVSTAT"]
            ),
            "CLNSIG": clinvar_rec.INFO["CLNSIG"],
            "gene_info": normalize_gene_info(clinvar_rec),
        }
        if "None" not in clinvar_dict["alt"] and clinvar_dict["id"] not in uniq_ids:
            # extract all unique clinvars that do not contain the None nucl as alternate
            clinvars_info.append(clinvar_dict)
            uniq_ids.append(clinvar_dict["id"])
    return clinvars_info


def convert_review_status2stars(
    clinvar_stars_df: pd.DataFrame,
    star_status2int: dict[str, int],
    clinvar_rev_status: list[str],
) -> int:
    """
    Convert CLNREVSTAT (review status) tab from clinvar vcf file to number of review stars
    for unknown description -> star= -1
    ClinVar review status documentation: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
    """
    rev_status = [review_elem.replace("_", " ") for review_elem in clinvar_rev_status]

    rev_status = ",".join(rev_status)
    if (
        rev_status not in clinvar_stars_df.Review_status.values
    ):  # if retrieved status not in status
        return star_status2int["unknown review status"]
    else:
        return star_status2int[
            clinvar_stars_df.loc[clinvar_stars_df["Review_status"] == rev_status][
                "Number_of_gold_stars"
            ].iloc[0]
        ]


def normalize_gene_info(clinvar_rec: dict) -> list:
    """
    Normalize ClinVar gene information
    Returns
    -------
    list of str, int
        list of genes containing [gene symbol, gene id]
    """

    logger.debug("Normalize gene information")
    if "GENEINFO" in clinvar_rec.INFO:
        genes_info = []
        for gene in clinvar_rec.INFO["GENEINFO"].strip().split("|"):
            genes_info.append(gene.split(":"))
        return genes_info
    else:
        return None


def uniq_new_filter(new_clinvars: list[dict[str, str]]) -> list[dict[str, str]]:
    """
    Filter clinvar records to keep the ones not already found in matched list
    """

    logger.debug("\nFilter all new clinvars to keep the unique ones only")
    ids, uniq_clinvars = [], []
    for clinvar in new_clinvars:
        if clinvar["id"] not in ids:
            uniq_clinvars.append(clinvar)
            ids.append(clinvar["id"])
    logger.debug(f"Unique clinvars: {uniq_clinvars}")
    return uniq_clinvars


def strand_filter(
    matched_clinvars: list[dict[str, str]], transcript_strand: str
) -> list[dict[str, str]]:
    """
    Filter matched ClinVars to ensure the same strand with affected transcript

    Returns
    -------
    list of dict of str: str
        filtered clinvars by strand information
    """

    logger.debug(f"\nFilter by strand the {len(matched_clinvars)} matching clinvars")
    filtered_clinvars = []
    for candidate_clinvar in matched_clinvars:
        if candidate_clinvar["gene_info"]:
            # get the symbol for the first registered gene in ClinVar
            gene_symbol, gene_id = candidate_clinvar["gene_info"][0]
        else:
            gene_symbol = None
        logger.debug(f"clinvar's gene: {gene_symbol}")

        ### ### ###
        # save a clinvar entry if
        # a) the strand is the same as for the PyEnsembl annotation
        # b) the significance field is filled
        ### ### ###
        if (
            get_clinvar_strand(hugo_genes_df, gene_symbol) == transcript_strand
            and "CLNSIG" in candidate_clinvar
        ):
            if "None" not in candidate_clinvar["alt"]:
                # do not extract a clinvar record that contains the None nucl as alternate
                filtered_clinvars.append(candidate_clinvar)
    logger.debug(f"Strand-filtered clinvars: {filtered_clinvars}")
    return filtered_clinvars


def get_clinvar_strand(hugo_genes_df: pd.DataFrame, gene_symbol: str) -> str:
    """
    Get strand for gene found in clinvar entry
    """

    logger.debug("Get strand for gene found in ClinVar entry")
    if gene_symbol in hugo_genes_df.index:
        return hugo_genes_df.loc[gene_symbol].strand
    else:
        return None


def cache_clinvars(
    chr: str, pos: int, strand: str, extracted_clinvars: list[dict[str, str]]
) -> None:
    """
    Save extracted clinvars to cache
    """

    logger.debug(f"Cache clinvars for pos: chr{chr},{pos},{strand}")
    if chr in clinvars_cache:
        if pos in clinvars_cache[chr]:
            if strand in clinvars_cache[chr][pos]:
                clinvars_cache[chr][pos][strand] = (
                    clinvars_cache[chr][pos][strand] + extracted_clinvars
                )
            else:
                clinvars_cache[chr][pos][strand] = extracted_clinvars
        else:
            clinvars_cache[chr][pos] = {strand: extracted_clinvars}
    else:
        clinvars_cache[chr] = {pos: {strand: extracted_clinvars}}


def construct_clinvar_prot_change(clinvar_rec, var_codon_info):
    """
    Construct protein change for clinvar record

    For termination codon the return protein change will contain the 3-letter protein letter 'Ter'
    as discussed: https://varnomen.hgvs.org/recommendations/protein/variant/frameshift/

    Parameters
    ----------
    clinvar_rec: dict of str : str
        clinvar record information
    var_codon_info: dict of str: str
        variant codon information

    Returns
    -------
    int
        protein change position
    str
        clinvar record's protein change
    """

    logger.debug(f"Construct protein change for clinvar record: {clinvar_rec}")
    # get position of clinvar nucleotide
    var_pos_idx = var_codon_info["genomic_pos"].index(clinvar_rec["pos"])
    logger.debug(f"ClinVar record is found in codon position: {var_pos_idx}")
    ### ### ###
    # construct the reference coding sequence for clinvar record
    ### ### ###
    clinvar_ref_seq, clinvar_alt_seq = [], []
    for idx, nucl in enumerate(var_codon_info["seq_ref"]):
        if idx == var_pos_idx:
            clinvar_ref_seq.append(clinvar_rec["ref"])
        else:
            clinvar_ref_seq.append(nucl)
    logger.debug(f"clinvar ref: {clinvar_ref_seq}")
    ### ### ###
    # and then construct the alternate sequence
    ### ### ###
    for idx, nucl in enumerate(var_codon_info["seq_ref"]):
        if idx == var_pos_idx:
            clinvar_alt_seq.append(clinvar_rec["alt"])
        else:
            clinvar_alt_seq.append(nucl)
    logger.debug(f"clinvar alt: {clinvar_alt_seq}")

    ### ### ###
    # translate these ref and alt sequence to protein edit
    ### ### ###
    clinvar_ref_translated = str(
        Seq(normalize_codon_nucleotides(clinvar_ref_seq)).translate()
    )
    clinvar_ref_aa = "".join(convert_1to3_aa(clinvar_ref_translated))
    clinvar_alt_translated = str(
        Seq(normalize_codon_nucleotides(clinvar_alt_seq)).translate()
    )
    clinvar_alt_aa = "".join(convert_1to3_aa(clinvar_alt_translated))
    logger.debug(
        f"ref_translated: {clinvar_ref_translated}, alt_translated: {clinvar_alt_translated}"
    )
    return (
        var_codon_info["prot_start"],
        clinvar_ref_aa + str(var_codon_info["prot_start"]) + clinvar_alt_aa,
    )


def compact_clinvar_entries(clinvar_records: list[dict]) -> str:
    """
    Compact ClinVar entries
    """

    logger.debug("Compact clinvar entries")
    clinvar_attributes_ordered = [
        "id",
        "pos",
        "ref",
        "alt",
        "CLNDISDB",
        "CLNREVSTAT",
        "CLNSIG",
    ]

    compacted_clinvars = []
    for clinvar in clinvar_records:
        clinvar_str = []
        for attr in clinvar_attributes_ordered:
            if attr != "CLNSIG":
                clinvar_str.append(attr + ":" + str(clinvar[attr]))
            else:
                clinvar_str.append(
                    attr
                    + ":"
                    + "::".join(
                        [str(signif).replace("_", " ") for signif in clinvar[attr]]
                    )
                )
        # save compacted clinvar and append to all compacts
        compacted_clinvar = ",".join(clinvar_str)
        logger.debug(f"compacted clinvar: {compacted_clinvar}")
        compacted_clinvars.append(compacted_clinvar)
    logger.debug(f"compacted clinvars: {compacted_clinvars}")
    return ";".join(compacted_clinvars)
