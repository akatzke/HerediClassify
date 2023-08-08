#!/usr/bin/env python3

import pyensembl
import logging

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
