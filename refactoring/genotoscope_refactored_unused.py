#!/usr/bin/env python3

import logging

import pyensembl
from variant import VariantInfo, TranscriptInfo

from genotoscope_exon_skipping import (
    get_transcript_exon_offsets,
    is_transcript_in_positive_strand,
    is_transcript_type_splice_acceptor_donor,
    parse_variant_intron_pos,
)

logger = logging.getLogger("GenOtoScope_Classify.tmp")


def find_exon_by_var_pos(
    ref_transcript: pyensembl.transcript.Transcript,
    transcript: TranscriptInfo,
    variant: VariantInfo,
    is_genomic: bool,
):
    """
    Find variant exon index by variant coding position, exon index is 1-based
    Returns
    -------
    list of int
        indices of exons overlapping with variant (1-based)
    int
        variant start position as offset in first overlapping exon
    int
        variant end position as offset in last overlapping exon
    """

    logger.debug("Find exon indices containing the variant")
    if not (transcript.exon or transcript.intron):
        logger.debug("Exon information is not found in VEP => use VEP spliced offset")
        is_genomic = False

    exon_positions = get_transcript_exon_offsets(ref_transcript, is_genomic)

    overlap_exon_indices = []
    var_start_offset, var_end_offset = (
        0,
        0,
    )  # save the variant start and end offset in the overlapping exon regions
    var_coding = transcript.var_hgvs
    ### ### ###
    # get the strand direction
    ### ### ###
    if is_transcript_in_positive_strand(ref_transcript):
        strand_direction = +1
    else:
        strand_direction = -1
    if is_genomic:
        if not is_transcript_type_splice_acceptor_donor(transcript.var_type):
            ### ### ###
            # Exonic variant
            ### ### ###
            logger.debug("Exonic variant type")
            var_start = variant.genomic_start
            var_end = variant.genomic_end
        else:
            logger.debug(
                "Intron variant type => update as variant pos the starting position of skipped exon"
            )
            split_symbols, intron_offsets, directions2exon = parse_variant_intron_pos(
                var_coding
            )

            ### ### ###
            # Intronic variant
            # Parse affected exon index by VEP (1-based)
            # add splice site direction to affected exon index
            # (currently modeled that) intron variant can disrupt at most one exon, so equal variant end with its start position
            ### ### ###

            if transcript.exon and transcript.intron:
                logger.debug(
                    "Transcript info contains both exon and intron offset => use intron offset"
                )
                exon_idx = transcript.intron - 1
            elif transcript.exon:
                exon_idx = transcript.exon - 1
            elif transcript.intron:
                exon_idx = transcript.intron - 1
            else:
                try:
                    assert transcript.exon or transcript.intron
                except AssertionError:
                    logger.error(
                        f"Variant does not contain exon or intron index in VEP column\n => variant position: {variant.to_string()}",
                        exc_info=True,
                    )
            if exon_idx + directions2exon[0] >= 0:
                skipped_exon = exon_positions[exon_idx + directions2exon[0]]
            else:
                logger.debug("Skipping first exon")
                skipped_exon = exon_positions[0]
                var_start = skipped_exon[0]
                var_end = var_start
                logger.debug(
                    f"Updated variant start:{var_start}, end: {var_end} on exon idx: {exon_idx + directions2exon[0]}"
                )
    else:
        # use cDNA offset for both frameshift and nonsense mutation
        var_start = int(var_coding.pos.start.base)
        var_end = int(var_coding.pos.end.base)
        if var_start < 0:
            var_start = 0
        if var_end < 0:
            if transcript.diff_len > 0:
                var_end = var_start + transcript.diff_len
            else:
                # deletion case
                var_end = var_start

        # VEP positions include only the coding sequence of the transcript,
        # so you need to add the length of the sequence of the first exons up to start codon (ATG)
        len_first_exons_up_start_codon = ref_transcript.start_codon_spliced_offsets[0]
        var_start = var_start + len_first_exons_up_start_codon
        var_end = var_end + len_first_exons_up_start_codon

    ### ### ###
    # find variant position into exon intervals
    ### ### ###
    for exon_idx, exon_interval in enumerate(exon_positions):
        if is_genomic:
            normalized_exon_interval = range(
                exon_interval[0],
                exon_interval[1] + 2 * strand_direction,
                strand_direction,
            )
        else:
            normalized_exon_interval = range(exon_interval[0], exon_interval[1] + 1)
            logger.debug(f"Exon interval: {normalized_exon_interval}")
            if var_start in normalized_exon_interval:
                overlap_exon_indices.append(exon_idx + 1)
                break
            logger.debug(
                f"Search var_start: {var_start} found in exon(s): {overlap_exon_indices}"
            )

    ### ### ###
    # find variant start and end offset in exons
    ### ### ###
    if len(overlap_exon_indices) == 1:  # variant included in only one exon
        if is_transcript_type_splice_acceptor_donor(transcript.var_type):
            # exon-skipping, get as offset the total range of skipped exon
            var_start_offset = 0
            # negative strand transcripts contain higher start position than end position
            var_end_offset = abs(
                exon_positions[overlap_exon_indices[0] - 1][1]
                - exon_positions[overlap_exon_indices[0] - 1][0]
            )
        else:
            if is_genomic:
                if strand_direction == 1:
                    # for positive strand, compute distance from exon lower position (start)
                    var_start_offset = (
                        var_start - exon_positions[overlap_exon_indices[0] - 1][0]
                    )
                    var_end_offset = (
                        var_end - exon_positions[overlap_exon_indices[0] - 1][0]
                    )
                else:
                    # for negative strand, compute distance from exon higher position (end)
                    var_start_offset = (
                        exon_positions[overlap_exon_indices[0] - 1][0] - var_end
                    )
                    var_end_offset = (
                        exon_positions[overlap_exon_indices[0] - 1][0] - var_start
                    )
            else:
                var_start_offset = (
                    var_start - exon_positions[overlap_exon_indices[0] - 1][0]
                )
                var_end_offset = (
                    var_end - exon_positions[overlap_exon_indices[0] - 1][0]
                )
    else:
        try:
            assert len(overlap_exon_indices) > 0
        except AssertionError:
            logger.error(
                f"Overlapping exon indices should be more than 0\n=> variant position: {variant.to_string()}",
                exc_info=True,
            )
        ### ### ###
        # variant included in more than one exons,
        # so start is in the first exon, end is in the last exon
        ### ### ###
        if is_genomic:
            if strand_direction == 1:
                var_start_offset = (
                    var_start - exon_positions[overlap_exon_indices[0] - 1][0]
                )
                var_end_offset = (
                    var_end
                    - exon_positions[
                        overlap_exon_indices[len(overlap_exon_indices)] - 1
                    ][0]
                )
            else:
                # for negative strand, use the highest value (exon start) to compute the offset for the start and exon position of the variant
                var_start_offset = (
                    exon_positions[overlap_exon_indices[0] - 1][0] - var_end
                )
                var_end_offset = (
                    exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][
                        0
                    ]
                    - var_start
                )
        else:
            var_start_offset = (
                var_start - exon_positions[overlap_exon_indices[0] - 1][0]
            )
            var_end_offset = (
                var_end
                - exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][0]
            )
    try:
        assert var_start_offset >= 0 and var_end_offset >= 0
    except AssertionError:
        logger.error(
            f"Variant start and end offset should be higher than 0\n=> variant position: {variant.to_string()}",
            exc_info=True,
        )
        logger.debug(
            f"Overlap exon indices: {overlap_exon_indices}, var start offset: {var_start_offset}, var end offset: {var_end_offset}"
        )
    return overlap_exon_indices, var_start_offset, var_end_offset
