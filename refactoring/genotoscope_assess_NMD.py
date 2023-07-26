#!/usr/bin/env python3

import logging
import pyensembl
from genotoscope_exon_skipping import (
    get_transcript_exon_offsets,
    is_transcript_in_positive_strand,
    assess_exon_skipping,
    find_exon_by_ref_pos,
)
from genotoscope_refactored_unused import find_exon_by_var_pos
from refactoring.variant import TranscriptInfo, VariantInfo

logger = logging.getLogger("GenOtoScope_Classify.PVS1.assess_NMD")


def assess_NMD(transcript: TranscriptInfo, variant: VariantInfo) -> tuple:
    """
    Examine if "non sense mediated mRNA decay" (NMD) will occur for current variant
    Following Tayoun et al. and Zhiyuan et al.  NMD is not predicted to occur if:
    a) premature termination codon (PTC) occurs in the last exon
    b) PTC occur in the (3') last 50 nucleotides of the penultimate exon
    c) transcript contains no introns
    d) PTC occurs inbetween the first 200 bases from start codon

    Returns
    -------
    bool
        NMD predicted to occur (True), not to occur (False)
    list of dict of str: str or int
        var containing exon genomic positions
    """

    logger.debug(f"Assess NMD for transript: {transcript.transcript_id}")
    NMD_occurs = True
    ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
        transcript.transcript_id
    )
    # construct exon positions using chromosome position
    exon_positions = get_transcript_exon_offsets(ref_transcript, True)
    num_exon_positions = len(exon_positions)
    ### ### ###
    # 4 cases: find exon containing the variant position
    # a) intronic variant that genomic start and/or end overlaps with exonic range
    # b) splice acceptor variant (that results to skipping)
    # c) exonic variant whose both start and stop overlap exonic range
    # d) exonic variant whose start or stop overlaps exonic range
    ### ### ###

    intron_variant_types = {
        "splice_acceptor_variant",
        "splice_acceptor",
        "splice_donor_variant",
        "splice_donor",
    }

    if transcript.var_type in intron_variant_types and (
        is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_start)
        or is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_end)
    ):
        ### ### ###
        # a) variant overlaps both intronic and exonic region
        ### ### ###
        logger.debug(
            "Variant genomic position overlaps both intronic and coding exonic range => use exon offsets (VEP coordinates) to find affected exon"
        )
        (
            exons_containing_var,
            var_exon_start_offset,
            var_exon_end_offset,
        ) = find_exon_by_var_pos(ref_transcript, transcript, variant, False)
        if is_transcript_type_splice_acceptor_donor(transcript.var_type):
            # if variant is splice acceptor, update exon skipping boolean flags
            (
                split_symbols,
                intron_offsets,
                directions2exon,
            ) = parse_variant_intron_pos(transcript.var_hgvs)
            (
                exon_containing_var,
                is_exon_skipped,
                var_exon_start_offset,
                var_exon_end_offset,
                variant_skips_start_codon_exon,
                variant_skips_stop_codon_exon,
                variant_skips_coding_exon,
            ) = assess_exon_skipping(
                ref_transcript, transcript, intron_offsets, variant
            )
    elif is_transcript_type_splice_acceptor_donor(transcript.var_type):
        logger.debug("assess NMD for splice acceptor variant")
        ### ### ###
        # b) splice acceptor
        # Currently handling only splice acceptor case that results to skipping
        # find affected exon index, start and stop and boolean flags by assess_exon_skipping()
        ### ### ###
        split_symbols, intron_offsets, directions2exon = parse_variant_intron_pos(
            transcript.var_hgvs
        )
        (
            exons_containing_var,
            is_exon_skipped,
            var_exon_start_offset,
            var_exon_end_offset,
            variant_skips_start_codon_exon,
            variant_skips_stop_codon_exon,
            variant_skips_coding_exon,
        ) = assess_exon_skipping(ref_transcript, transcript, intron_offsets, variant)
        logger.debug("")
    elif not is_transcript_type_splice_acceptor_donor(transcript.var_type) and (
        is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_start)
        and is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_end)
    ):
        ### ### ###
        # c) exonic variant with both genomic start and stop in exonic range
        ### ### ###
        logger.debug(
            "Exonic type of variant and its genomic start and stop both overlap exonic range"
        )
        (
            exons_containing_var,
            var_exon_start_offset,
            var_exon_end_offset,
        ) = find_exon_by_var_pos(ref_transcript, transcript, variant, True)
    elif not is_transcript_type_splice_acceptor_donor(transcript.var_type) and (
        is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_start)
        or is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_end)
    ):
        ### ### ###
        # d) exonic variant with either genomic start or stop in exonic range
        ### ### ###
        logger.debug(
            "Exonic type of variant and its genomic start or stop overlaps exonic range"
        )
        (
            exons_containing_var,
            var_exon_start_offset,
            var_exon_end_offset,
        ) = find_exon_by_var_pos(ref_transcript, transcript, variant, False)

    ### ### ### ### ### ###
    # find interacting exons genomic positions
    # if splice acceptor variant results to exon skipping
    # or the variant is exonic
    ### ### ### ### ### ###

    exon_variant_types = {
        "stop_gained",
        "stop_lost",
        "frameshift",
        "frameshift_variant",
        "inframe_insertion",
        "inframe_deletion",
        "disruptive_inframe_insertion",
        "disruptive_inframe_deletion",
        "conservative_inframe_insertion",
        "conservative_inframe_deletion",
    }

    if transcript.coding_exon_skipped or transcript.var_type in exon_variant_types:
        logger.debug("Variant applicable for stop codon searching")
        if (
            transcript.start_codon_exon_skipped and transcript.stop_codon_exon_skipped
        ) or not transcript.start_codon_exon_skipped:
            logger.debug(
                "Update exon positions for variant that skips both start and stop exons or does not skip start codon exon"
            )
            affected_exons_pos = find_affected_exons_pos(
                ref_transcript,
                exons_containing_var,
                var_exon_start_offset,
                var_exon_end_offset,
                variant,
            )

            exon2stop_index, num_exons = search_stop_codon(
                ref_transcript,
                transcript,
                exons_containing_var,
                is_exon_skipped,
                variant_skips_stop_codon_exon,
                variant,
            )

            if num_exons == 1:
                # intronless transcript
                NMD_occurs, NMD_comment = False, "Single exon"
            else:
                exons_with_stop_codon = sorted(list(exon2stop_index.keys()))
                if len(exon2stop_index) == 0:
                    # PTC not found in any exon
                    # => presumably in the 3' non coding region of the last exon
                    NMD_occurs, NMD_comment = False, "PTC after reference stop codon"
                elif (
                    exons_with_stop_codon[0] == num_exon_positions - 1
                    or exons_with_stop_codon[0] == num_exon_positions
                ):
                    # PTC in 3'-most 50 bases of penultimate or ultimate exon
                    NMD_occurs, NMD_comment = (
                        False,
                        "PTC in last exon or 3'-most 50 bases of penultimate exon",
                    )
                elif exon2stop_index[exons_with_stop_codon[0]] < 200:
                    # PTC in the first 200 bases from start codon
                    NMD_occurs, NMD_comment = (
                        False,
                        "PTC distance from start codon < 200",
                    )
                elif exon2stop_index[exons_with_stop_codon[0]] > 200:
                    # PTC after the first 200 bases from start codon
                    # but before the last exon junction complex EJC
                    NMD_occurs, NMD_comment = True, "PTC before last EJC"

            logger.debug(f"Transcript contains in total: {num_exons} exon(s)")
            logger.debug(f"positions of stop codons: {exon2stop_index}")
            logger.debug(
                f"NMD is predicted to occur: {NMD_occurs}, comment: {NMD_comment}"
            )
        else:
            affected_exons_pos = []
            NMD_occurs = False
    else:
        logger.debug("Variant type not applicable for stop codon search")
        affected_exons_pos = []

    return (
        NMD_occurs,
        affected_exons_pos,
    )


def is_genomic_pos_in_coding_exon(
    transcript: pyensembl.transcript.Transcript, genomic_pos: int
) -> bool:
    """
    Examine if genomic position is contained in coding exon
    """

    pos_in_coding_exon = False
    logger.debug(
        f"Examine if genomic pos: {genomic_pos} is contained in coding exon sequence"
    )
    if is_transcript_in_positive_strand(transcript):
        strand_direction = +1
    else:
        strand_direction = -1
    num_coding_exons = len(transcript.coding_sequence_position_ranges)

    ### ### ###
    # loop through values and check if variant overlap a coding exonic range
    ### ### ###
    for coding_exon_idx, coding_interval in enumerate(
        transcript.coding_sequence_position_ranges
    ):
        # set up exon coding start and end positions
        if strand_direction == 1:
            exon_coding_start = coding_interval[0]
            exon_coding_end = coding_interval[1]
        else:
            exon_coding_start = coding_interval[1]
            if coding_exon_idx + 1 == num_coding_exons:
                # update end of last exon to be the lowest chromosomal position of the stop codon
                exon_coding_end = transcript.stop_codon_positions[0]
            else:
                exon_coding_end = coding_interval[0]
        normalized_coding_interval = range(
            exon_coding_start, exon_coding_end + 2 * strand_direction, strand_direction
        )

        logger.debug(f"normalized interval: {normalized_coding_interval}")
        if genomic_pos in normalized_coding_interval:
            logger.debug(f"position in exon with offset: {coding_exon_idx}")
            pos_in_coding_exon = True
            break
    return pos_in_coding_exon


def search_stop_codon(
    ref_transcript: pyensembl.transcript.Transcript,
    transcript: TranscriptInfo,
    exons_containing_var,
    variant: VariantInfo,
):
    """
    Search for stop codon on all observed exonic coding sequences

    Returns
    -------
    dict of int: int
        map of exon index to (left-most) termination codon position
    int
        number of exons
    """

    logger.debug("Search stop codon over all exons' observed coding sequence")
    # create the variant coding sequence of the trancript

    var_coding_seq, diff_len = (
        transcript.var_seq,
        transcript.diff_len,
    )
    logger.debug("Difference of observed to reference, in length: {}".format(diff_len))
    current_codon_length, remain_codon_length, last_start = 0, 0, 0
    if is_transcript_in_positive_strand(ref_transcript):
        start_codon_first_pos = ref_transcript.start_codon_positions[0]
        # find index of exon that contains stop codon
        logger.debug(
            "Stop codon positions: {}".format(ref_transcript.stop_codon_positions)
        )
        stop_codon_first_base_exon_idx, _ = find_exon_by_ref_pos(
            ref_transcript, ref_transcript.stop_codon_positions[0], True
        )
        stop_codon_last_base_exon_idx, _ = find_exon_by_ref_pos(
            ref_transcript, ref_transcript.stop_codon_positions[2], True
        )
    else:
        start_codon_first_pos = ref_transcript.start_codon_positions[-1]
        # find index of exon that contains stop codon
        stop_codon_first_base_exon_idx, _ = find_exon_by_ref_pos(
            transcript, ref_transcript.stop_codon_positions[-1], True
        )
        stop_codon_last_base_exon_idx, _ = find_exon_by_ref_pos(
            transcript, ref_transcript.stop_codon_positions[0], True
        )

    # if stop codon starts on one exon and finishes on the next one
    # keep the last exon as the position of the stop codon
    stop_codon_exon_idx = max(
        stop_codon_first_base_exon_idx, stop_codon_last_base_exon_idx
    )
    exon_contains_coding_seq, exon_contains_start_codon, exon_contains_stop_codon = (
        False,
        False,
        False,
    )
    sum_obs_exon_coding_seq = 0
    exon2termination_codon_pos = {}
    exon_idx = 0
    num_exons = len(ref_transcript.exon_intervals)

    while not exon_contains_stop_codon and exon_idx < num_exons:
        (exon_start, exon_end) = ref_transcript.exon_intervals[exon_idx]
        logger.debug(
            "Exon idx: {}, start={} end={}".format(exon_idx + 1, exon_start, exon_end)
        )

        # examine if exon contains stop codon
        if exon_idx == stop_codon_exon_idx:
            exon_contains_stop_codon = True
        else:
            exon_contains_stop_codon = False

        ### ### ### ### ### ###
        # skip current exon if intron variant results to skip exon
        ### ### ### ### ### ###
        if exon_idx + 1 == exons_containing_var[0] and transcript.are_exons_skipped:
            logger.debug("Exon with idx:{} is skipped".format(exon_idx + 1))
            # we should not find that exon with start codon is skipped
            # because this case is captured by start_lost variant type
            exon_idx += 1
            continue

        ### ### ### ### ### ###
        # do not construct coding sequence for exons up to the exon that contains start codon
        ### ### ### ### ### ###
        if not exon_contains_coding_seq:
            if exon_start <= start_codon_first_pos <= exon_end:
                logger.debug(
                    "First exon with start codon has index: {}".format(exon_idx + 1)
                )
                logger.debug("last start: {}".format(last_start))
                exon_contains_start_codon = True
                exon_contains_coding_seq = True

                (
                    exon_contains_coding_seq,
                    exon_contains_start_codon,
                    last_start,
                    current_codon_length,
                    remain_codon_length,
                ) = construct_observed_exon_seq(
                    transcript,
                    exon_idx,
                    exons_containing_var,
                    exon_contains_start_codon,
                    exon_contains_stop_codon,
                    exon_contains_coding_seq,
                    diff_len,
                    last_start,
                    current_codon_length,
                    remain_codon_length,
                )
            else:
                logger.debug("UTR region in exon, before start codon")
                exon_idx += 1
                continue
        else:
            # exon contains coding sequence
            # construct current exon observed coding sequence
            (
                exon_contains_coding_seq,
                exon_contains_start_codon,
                last_start,
                current_codon_length,
                remain_codon_length,
            ) = construct_observed_exon_seq(
                transcript,
                exon_idx,
                exons_containing_var,
                exon_contains_start_codon,
                exon_contains_stop_codon,
                exon_contains_coding_seq,
                diff_len,
                last_start,
                current_codon_length,
                remain_codon_length,
            )

        if exon_contains_coding_seq:
            ### ### ### ### ### ###
            # get observed current exon coding sequence and codons
            ### ### ### ### ### ###
            if exon_contains_stop_codon and (
                exon_start <= start_codon_first_pos <= exon_end
            ):
                # add all remaining sequence if current exon contains both start and stop codons
                logger.debug("Exon contains both start and stop codons")
                obs_exon_coding_seq = var_coding_seq[last_start : len(var_coding_seq)]
                logger.debug(
                    "observed exonic codons seq= {}, length of seq= {}".format(
                        obs_exon_coding_seq, len(obs_exon_coding_seq)
                    )
                )
            # self.logger.debug("remaining: {}".format(
            # 	var_coding_seq[last_start + current_codon_length + remain_codon_length:len(var_coding_seq)]))
            elif exon_contains_stop_codon:
                # on exon containing the stop codon, add the remaining coding sequence to the observed codon sequence
                logger.debug("Exon contains stop codon")
                obs_exon_coding_seq = var_coding_seq[
                    last_start : last_start + current_codon_length + remain_codon_length
                ]
                logger.debug(
                    "observed exonic codons seq= {}, length of seq= {}".format(
                        obs_exon_coding_seq, len(obs_exon_coding_seq)
                    )
                )
                logger.debug(
                    "remaining_ending: {}".format(
                        var_coding_seq[
                            last_start
                            + current_codon_length
                            + remain_codon_length : len(var_coding_seq)
                        ]
                    )
                )
                logger.debug("transcript type: {}".format(transcript.var_type))
            # obs_exon_coding_seq = obs_exon_coding_seq + var_coding_seq[last_start+current_codon_length+remain_codon_length:len(var_coding_seq)]
            else:
                # update observed coding sequence for coding exon
                obs_exon_coding_seq = var_coding_seq[
                    last_start : last_start + current_codon_length
                ]
                logger.debug(
                    "observed exonic codons seq= {}, length of seq= {}".format(
                        obs_exon_coding_seq, len(obs_exon_coding_seq)
                    )
                )
            # self.logger.debug("remaining: {}".format(
            # 	var_coding_seq[last_start + current_codon_length:len(var_coding_seq)]))

            observed_codons = extract_codons(obs_exon_coding_seq)

            if exon_idx == stop_codon_exon_idx - 1:
                logger.debug("Process penultimate exon")
                # if exon_contains_stop_codon:
                # search termination codon on penultimate exon
                (
                    termination_codon_exists,
                    termination_codon_index,
                ) = search_termination_codon(observed_codons, True)
            else:
                # search termination codon on any other exon
                (
                    termination_codon_exists,
                    termination_codon_index,
                ) = search_termination_codon(observed_codons, False)

            if termination_codon_exists:
                # to create absolute termination codon indices
                # add the current sum of observed coding sequence
                exon2termination_codon_pos[exon_idx + 1] = (
                    sum_obs_exon_coding_seq + termination_codon_index * 3
                )
            sum_obs_exon_coding_seq += len(obs_exon_coding_seq)
        # update exon index
        exon_idx += 1

    logger.debug("sum observed exon coding seq: {}".format(sum_obs_exon_coding_seq))
    logger.debug(
        "length of var_coding_seq={} and sum processed coding seq={}".format(
            len(var_coding_seq), sum_obs_exon_coding_seq
        )
    )

    if not transcript.stop_codon_exon_skipped:
        try:
            assert sum_obs_exon_coding_seq == len(
                var_coding_seq
            ) or sum_obs_exon_coding_seq + len(var_coding_seq) - (
                last_start + current_codon_length + remain_codon_length
            ) == len(
                var_coding_seq
            )
        except AssertionError:
            logger.error(
                "Sum of observed exonic codons should equal the total length of the observed coding sequence\n=> variant position: {}".format(
                    variant.to_string()
                ),
                exc_info=True,
            )
    return exon2termination_codon_pos, num_exons
