#!/usr/bin/env python3

import pyensembl
import logging
import pandas as pd

from refactoring.variant import VariantInfo
from refactoring.genotoscope_exon_skipping import (
    is_transcript_in_positive_strand,
)

logger = logging.getLogger("GenOtoScope_Classify.unused_functions")


def convert_genomic2coding_pos(
    ref_transcript: pyensembl.transcript.Transcript,
    genomic_pos: int,
    variant: VariantInfo,
) -> int:
    """
    Convert genomic position into coding position
    """

    logger.debug("Convert variant chromosomal position to coding position")

    if is_transcript_in_positive_strand(ref_transcript):
        strand_direction = +1
        start_codon_first_pos = ref_transcript.start_codon_positions[0]
    else:
        strand_direction = -1
        start_codon_first_pos = ref_transcript.start_codon_positions[-1]

    sum_coding_len, pos_exon_offset = 0, -1
    logger.debug(f"start codon: {ref_transcript.start_codon_positions}")
    logger.debug(f"stop codon: {ref_transcript.stop_codon_positions}")
    num_coding_exons = len(ref_transcript.coding_sequence_position_ranges)

    ### ### ###
    # loop through coding sequence of exons
    # if variant chromosomal position is not included on current exon, add up of coding length
    # else compute the offset of the variant position from the current exon start (respecting strand direction)
    ### ### ###
    for coding_exon_idx, coding_interval in enumerate(
        ref_transcript.coding_sequence_position_ranges
    ):
        # set up exon coding start and end positions
        if strand_direction == 1:
            exon_coding_start = coding_interval[0]
            exon_coding_end = coding_interval[1]
        else:
            exon_coding_start = coding_interval[1]
            if coding_exon_idx + 1 == num_coding_exons:
                # update end of last exon to be the lowest chromosomal position of the stop codon
                exon_coding_end = ref_transcript.stop_codon_positions[0]
            else:
                exon_coding_end = coding_interval[0]
        normalized_coding_interval = range(
            exon_coding_start, exon_coding_end + strand_direction, strand_direction
        )
        logger.debug(f"normalized interval: {normalized_coding_interval}")
        if genomic_pos in normalized_coding_interval:
            logger.debug("position in exon")
            # find offset of genomic position in containing exon
            if strand_direction == +1:
                pos_exon_offset = genomic_pos - exon_coding_start
            else:
                pos_exon_offset = exon_coding_start - genomic_pos
            logger.debug(f"exon offset: {pos_exon_offset}")
            break
        else:
            # add coding length for each coding exon before exon containing genomic position
            sum_coding_len = (
                sum_coding_len + abs(exon_coding_end - exon_coding_start) + 1
            )
            logger.debug(f"sum of coding: {sum_coding_len}")

    if pos_exon_offset == -1:
        logger.debug("Normalize position that is before the start codon")
        # check if variant contains non-translated positions before the start codon
        # update position to be the first coding position
        sum_coding_len = 0
        if strand_direction == +1:
            if genomic_pos < start_codon_first_pos:
                pos_exon_offset = 0
        else:
            if genomic_pos > start_codon_first_pos:
                pos_exon_offset = 0
    try:
        assert pos_exon_offset != -1
    except AssertionError:
        logger.error(
            f"Chromosome position should be on any coding exon\n => variant position: {variant.to_string()}",
            exc_info=True,
        )
    return sum_coding_len + pos_exon_offset


def extract_pathogenic_var_upstream_closest_start_codon(
    ref_transcript, closest_start_codon_index, variant
):
    """
    Parameters
    ----------
    ref_transcript : pyensembl.transcript
        transcript object that variant lies in
    closest_start_codon_index : int
        closest start codon index
    variant : VariantInfo
        variant basic info

    Returns
    -------
    list of dict of str : int or str or list of str
        pathogenic clinvars between affected start codon and closest in-frame start codon
    """

    # print("Extract pathogenic variants in upstream of the closest codon")
    is_genomic = False
    # find genomic position for the original start codon
    (
        original_start_codon_exon_index,
        original_start_codon_exon_offset,
    ) = find_exon_by_ref_pos(ref_transcript, 1, is_genomic)

    # based on transcript directionality, add start,end co-ordinates
    if is_transcript_in_positive_strand(ref_transcript):
        transcript_strand = "+"
    else:
        transcript_strand = "-"

    original_start_codon_start = ref_transcript.start_codon_positions[0]
    original_start_codon_end = ref_transcript.start_codon_positions[-1]

    original_start_codon = {
        "exon_idx": original_start_codon_exon_index,
        "start": original_start_codon_start,
        "end": original_start_codon_end,
    }
    # print("original start codon: {}".format(original_start_codon))

    # find exon position of the closest in-frame start codon
    closest_start_codon_pos = (closest_start_codon_index + 1) * 3
    (
        closest_start_codon_exon_index,
        closest_start_codon_exon_offset,
    ) = find_exon_by_ref_pos(ref_transcript, closest_start_codon_pos, is_genomic)
    # convert exon position in genomic position for variants retrieval
    closest_start_codon_pos = convert_exon_pos2genomic_pos(
        ref_transcript, closest_start_codon_exon_index, closest_start_codon_exon_offset
    )
    # print("closet start codon pos: {}".format(closest_start_codon_pos))

    # based on transcript directionality add start,end co-ordinates
    if transcript_strand == "+":
        closest_start_codon_start = closest_start_codon_pos
        closest_start_codon_end = closest_start_codon_pos + 2
    else:
        # closest_start_codon_start = closest_start_codon_pos
        # closest_start_codon_end = transcript.exons[closest_start_codon_exon_index].to_dict()["end"]
        closest_start_codon_start = closest_start_codon_pos - 2
        closest_start_codon_end = closest_start_codon_pos

    closest_start_codon = {
        "exon_idx": closest_start_codon_exon_index,
        "start": closest_start_codon_start,
        "end": closest_start_codon_end,
    }
    # print("closest start codon: {}".format(closest_start_codon))

    # assert the start and end positions in the variant search ranges
    if transcript_strand == "+":
        try:
            assert original_start_codon["start"] < closest_start_codon["end"]
        except AssertionError:
            logger.error(
                f"Original and alternative start on the same exon; in positive strand transcript, original's start position should be lower than closest's end position\n=> variant position: {variant.to_string()}",
                exc_info=True,
            )
    else:
        try:
            assert closest_start_codon["start"] < original_start_codon["end"]
        except AssertionError:
            logger.error(
                f"Original and alternative start on the same exon; in negative strand transcript, original's start position should be higher than closest's end position\n=> variant position: {variant.to_string()}",
                exc_info=True,
            )
    start_codon2closest_start_ranges = set_up_start_end_pos(
        ref_transcript, original_start_codon, closest_start_codon
    )
    logger.debug("start codon ranges: {}".format(start_codon2closest_start_ranges))
    clinvars = extract_clinvars(
        start_codon2closest_start_ranges,
        ref_transcript.exons[0].to_dict()["strand"],
        variant,
    )
    if clinvars:
        filtered_clinvars = quality_filter_clinvars(clinvars, min_review_stars)
        return extract_pathogenic_clinvars(filtered_clinvars)
    else:
        # if there no extracted clinvars, we can't filter neither select the pathogenic ones
        # thus return None as pathogenic clinvar entries
        return None


def set_up_start_end_pos(
    ref_transcript, original_start_codon, closest_start_codon
) -> list[list[int]]:
    """
    Set up start and end positions that are between the start codon and the closest in-frame start codon

    Parameters
    ----------
    ref_transcript : pyensembl.transcript
        transcript object that positions lie in
    original_start_codon : dict of str: str or int
        original start codon information
    closest_start_codon:
        closest inframe start codon information

    Returns
    -------
    list of list of int
        start and end genomic position for all exons between the start codon and closest in-frame start codon
    """

    start_end = []
    if original_start_codon["exon_idx"] == closest_start_codon["exon_idx"]:
        # original start and closest in-frame start codon are on the same exon
        if is_transcript_in_positive_strand(ref_transcript):
            start = original_start_codon["start"]
            end = closest_start_codon["end"]
        else:
            start = closest_start_codon["start"]
            end = original_start_codon["end"]
        start_end.append([start, end])
    elif closest_start_codon["exon_idx"] > original_start_codon["exon_idx"]:
        # closest start codon lies on exon on the downstream of the exon containing the original start codon
        if is_transcript_in_positive_strand(ref_transcript):
            transcript_strand = "+"
            start_end.append(
                [
                    original_start_codon["start"],
                    ref_transcript.exons[original_start_codon["exon_idx"]].to_dict()[
                        "end"
                    ],
                ]
            )
        else:
            transcript_strand = "-"
            start_end.append(
                [
                    ref_transcript.exons[original_start_codon["exon_idx"]].to_dict()[
                        "start"
                    ],
                    original_start_codon["end"],
                ]
            )
        for overlap_exon_idx in range(
            original_start_codon["exon_idx"] + 1, closest_start_codon["exon_idx"]
        ):
            overlap_exon = ref_transcript.exons[overlap_exon_idx]
            start_end.append(
                [overlap_exon.to_dict()["start"], overlap_exon.to_dict()["end"]]
            )

        if transcript_strand == "+":
            start_end.append(
                [
                    ref_transcript.exons[closest_start_codon["exon_idx"]].to_dict()[
                        "start"
                    ],
                    closest_start_codon["start"],
                ]
            )
        else:
            start_end.append(
                [
                    closest_start_codon["end"],
                    ref_transcript.exons[closest_start_codon["exon_idx"]].to_dict()[
                        "end"
                    ],
                ]
            )
    else:
        try:
            assert 1 == 0
        except AssertionError:
            logger.error(
                f"AssertError: Closest start codon lies in an exon before the original start codon\n=> original start codon: {original_start_codon}",
                exc_info=True,
            )

    for positions in start_end:
        assert positions[0] <= positions[1], "Start should be lower than end position"
    return start_end


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
