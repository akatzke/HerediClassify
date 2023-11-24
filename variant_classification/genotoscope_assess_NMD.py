#!/usr/bin/env python3

import logging
import pyensembl

from genotoscope_exon_skipping import (
    get_transcript_exon_offsets,
    is_transcript_in_positive_strand,
    find_exon_by_ref_pos,
    find_exon_by_var_pos,
)
from variant import TranscriptInfo, VariantInfo
from var_type import VARTYPE_GROUPS


logger = logging.getLogger("GenOtoScope_Classify.PVS1.assess_NMD")


def assess_NMD_threshold(
    transcript: TranscriptInfo, threshold: int, clin_transcript: str
) -> tuple:
    """
    Examine if position of variant is located before or after given threshold for NMD
    """
    if transcript.transcript_id == clin_transcript:
        if transcript.var_start >= threshold:
            NMD_affected_exon = []
            return True, NMD_affected_exon
        else:
            return False, []
    else:
        return True, []


def assess_NMD_intronic_variant(
    transcript: TranscriptInfo,
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    are_exons_skipped: bool,
    skipped_exons: list[int],
    start_codon_exon_skipped: bool,
    stop_codon_exon_skipped,
    coding_exon_skipped: bool,
    var_coding_seq: str,
    diff_len: int,
) -> tuple:
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

    # construct exon positions using chromosome position
    exon_positions = get_transcript_exon_offsets(ref_transcript, True)
    num_exon_positions = len(exon_positions)

    (skipped_exons, var_exon_start_offset, var_exon_end_offset) = find_exon_by_var_pos(
        ref_transcript, transcript, variant, False, diff_len
    )
    NMD_comment = "NMD not predicted"
    if coding_exon_skipped:
        logger.debug("Variant applicable for stop codon searching")
        if (
            start_codon_exon_skipped and stop_codon_exon_skipped
        ) or not start_codon_exon_skipped:
            logger.debug(
                "Update exon positions for variant that skips both start and stop exons or does not skip start codon exon"
            )
            affected_exons_pos = find_pos_affected_exons(
                ref_transcript,
                skipped_exons,
                variant,
            )

            exon2stop_index, num_exons = search_stop_codon(
                transcript,
                variant,
                ref_transcript,
                are_exons_skipped,
                skipped_exons,
                stop_codon_exon_skipped,
                var_coding_seq,
                diff_len,
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


def assess_NMD_exonic_variant(
    transcript: TranscriptInfo,
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    var_coding_seq: str,
    diff_len: int,
) -> tuple:
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

    # construct exon positions using chromosome position
    exon_positions = get_transcript_exon_offsets(ref_transcript, True)
    num_exon_positions = len(exon_positions)
    if is_genomic_pos_in_coding_exon(
        ref_transcript, variant.genomic_start
    ) and is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_end):
        # exonic variant with both genomic start and stop in exonic range
        logger.debug(
            "Exonic type of variant and its genomic start and stop both overlap exonic range"
        )
        (
            exons_containing_var,
            var_exon_start_offset,
            var_exon_end_offset,
        ) = find_exon_by_var_pos(ref_transcript, transcript, variant, True, diff_len)
    elif is_genomic_pos_in_coding_exon(
        ref_transcript, variant.genomic_start
    ) or is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_end):
        # exonic variant with either genomic start or stop in exonic range
        logger.debug(
            "Exonic type of variant and its genomic start or stop overlaps exonic range"
        )
        (
            exons_containing_var,
            var_exon_start_offset,
            var_exon_end_offset,
        ) = find_exon_by_var_pos(ref_transcript, transcript, variant, False, diff_len)
    else:
        raise ValueError(
            "Genomic position in transcript coding region expected, but variant no located in coding region."
        )

    NMD_comment = "NMD is not predicted"

    if any(
        var_type_exonic.value == var_type.value
        for var_type_exonic in VARTYPE_GROUPS.EXONIC.value
        for var_type in transcript.var_type
    ):
        logger.debug("Variant applicable for stop codon searching")
        logger.debug(
            "Update exon positions for variant that skips both start and stop exons or does not skip start codon exon"
        )
        affected_exons_pos = find_pos_affected_exons(
            ref_transcript,
            exons_containing_var,
            variant,
        )

        exon2stop_index, num_exons = search_stop_codon(
            transcript,
            variant,
            ref_transcript,
            False,
            exons_containing_var,
            False,
            var_coding_seq,
            diff_len,
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
        logger.debug(f"NMD is predicted to occur: {NMD_occurs}, comment: {NMD_comment}")
    else:
        logger.debug("Variant type not applicable for stop codon search")
        NMD_occurs = False
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


def find_pos_affected_exons(
    ref_transcript: pyensembl.transcript.Transcript,
    affected_exon_idxs: list[int],
    variant: VariantInfo,
) -> list[dict]:
    """
    Find start and end position of affected exons
    """

    logger.debug("Find affected exon positions")

    interact_exon_pos = []
    for affected_exon in affected_exon_idxs:
        ref_exon = ref_transcript.exons[affected_exon - 1]
        logger.debug(f"selected exon: {ref_exon.id}")
        interact_exon_pos.append(
            {
                "exon_id": ref_exon.id,
                "exon_start": ref_exon.start,
                "exon_end": ref_exon.end,
            }
        )
    # assert that created information for interacting exon have start <= end
    for interacting_exon in interact_exon_pos:
        try:
            assert interacting_exon["exon_start"] <= interacting_exon["exon_end"]
        except AssertionError:
            logger.error(
                f"Affected exon start is higher than end\n=> variant position: {variant.to_string}",
                exc_info=True,
            )

    logger.debug(f"Genomic positions of affected exons: {interact_exon_pos}")
    return interact_exon_pos


def search_stop_codon(
    transcript: TranscriptInfo,
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    are_exons_skipped: bool,
    exons_containing_var: list[int],
    stop_codon_exon_skipped: bool,
    var_coding_seq: str,
    diff_len: int,
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
            ref_transcript, ref_transcript.stop_codon_positions[-1], True
        )
        stop_codon_last_base_exon_idx, _ = find_exon_by_ref_pos(
            ref_transcript, ref_transcript.stop_codon_positions[0], True
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
        if exon_idx + 1 == exons_containing_var[0] and are_exons_skipped:
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
                    ref_transcript,
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
                ref_transcript,
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

    if not stop_codon_exon_skipped:
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


def search_termination_codon(
    exon_codons: list[str], is_penultimate_exon: bool
) -> tuple:
    """
    Search termination codon in exon codons

    Returns
    -------
    bool
        termination codon exists in exon codons (True), otherwise it doesn't (False)
    int
        termination codon index (1-index), if termination codon doesn't exist equals the length of the search coding region
    """

    if is_penultimate_exon:
        logger.debug("Search termination codon in penultimate exon")
    termination_codon_exists = False
    if is_penultimate_exon:
        # start from the most 3' 50 base ~= 50/3 = 17 codons
        search_start = len(exon_codons) - 17
        sequence_limit = "the most 3' 50 bases = 17 codons"
    # print("Current exon observed sequence ")
    else:
        # not penultimate exon
        search_start = 0
        sequence_limit = "all sequence"
    logger.debug(
        f"Current exon observed sequence ({sequence_limit}) is: {exon_codons[search_start : len(exon_codons)]}"
    )
    # find the left most (5'-most) termination codon in sequence's codons
    term_codons_indices = [
        get_codon_index(exon_codons[search_start : len(exon_codons)], term_codon)
        for term_codon in ["TAG", "TAA", "TGA"]
    ]
    left_most_term_index = min(term_codons_indices)

    # update termination codon existence flag and position
    if left_most_term_index == len(exon_codons[search_start : len(exon_codons)]):
        termination_codon_exists = False
        termination_codon_index = -1
    else:
        termination_codon_exists = True
        termination_codon_index = left_most_term_index
    logger.debug(f"Termination codon found: {termination_codon_exists}")
    return termination_codon_exists, termination_codon_index


def get_codon_index(seq_codons: list[str], target_codon: str) -> int:
    """
    Search target codon on sequence codons
    and return its index
    """
    try:
        codon_index = seq_codons.index(target_codon) + 1
    except ValueError:
        codon_index = len(seq_codons)
    return codon_index


def extract_codons(sequence: str) -> list[str]:
    """
    Extract codons from given sequence, which start where the frame start offset is
    """
    codons = []
    for pos in range(0, len(sequence), 3):
        codons.append(sequence[pos : pos + 3])
    return codons


def construct_observed_exon_seq(
    ref_transcript: pyensembl.transcript.Transcript,
    exon_idx: int,
    exons_containing_var: list[int],
    exon_contains_start_codon: bool,
    exon_contains_stop_codon: bool,
    exon_contains_coding_seq: bool,
    diff_len: int,
    last_start: int,
    current_codon_length: int,
    remain_codon_length: int,
):
    """
    Construct observed exon coding sequence

    Returns
    -------
    exon_contains_start_codon : bool
        exon contains start codon (True), otherwise exon contains only UTR (False)
    exon_contains_coding_seq : bool
        exon contains coding sequence (True), otherwise (False)
    last_start : int
        updated last coding exonic start position
    current_codon_length : int
        updated current exon coding length
    remain_codon_length : int
        updated coding length remained from previous (left) exon
    """

    logger.debug("Construct observed exon sequence")
    logger.debug(f"Diff len: {diff_len}")
    exon = ref_transcript.exons[exon_idx]
    exon_start, exon_end = exon.to_dict()["start"], exon.to_dict()["end"]
    # print("exon start: {}, exon end: {}".format(exon_start, exon_end))

    # set up first position of start codon based on transcript strand
    if is_transcript_in_positive_strand(ref_transcript):
        start_codon_first_pos = ref_transcript.start_codon_positions[0]
    else:
        start_codon_first_pos = ref_transcript.start_codon_positions[-1]

    # set up stop codon as the end of the last exon
    if exon_contains_stop_codon:
        logger.debug("this is last exon, update its end to be the stop codon position")
        logger.debug("before: exon start:{}, exon end:{}".format(exon_start, exon_end))
        # exon_end = transcript.stop_codon_positions[-1]
        if is_transcript_in_positive_strand(ref_transcript):
            exon_start = exon_start
            exon_end = ref_transcript.stop_codon_positions[-1]
        else:
            # print("stop codon positions: {}".format(transcript.stop_codon_positions))
            exon_start = ref_transcript.stop_codon_positions[0]
            exon_end = exon_end
    logger.debug(f"updated: exon start:{exon_start}, exon end:{exon_end}")
    coding_exonic_length = 0

    ### ### ### ### ### ###
    # update current exonic coding sequence
    # and its remaining for the next exon
    ### ### ### ### ### ###
    if exon_contains_coding_seq and exon_contains_start_codon:
        # first exon with coding sequence
        logger.debug("First exon containing coding sequence")
        logger.debug(f"last start: {last_start}".format(last_start))
        if is_transcript_in_positive_strand(ref_transcript):
            coding_exonic_length = (exon_end - start_codon_first_pos) + 1
        else:
            coding_exonic_length = (start_codon_first_pos - exon_start) + 1
        logger.debug(
            f"before: adding variant edit, coding_exonic_length={coding_exonic_length}"
        )
        if exon_idx + 1 == exons_containing_var and diff_len != 0:
            logger.debug(
                f"add variant edit len={diff_len} on the current length of exon idx:{exon_idx + 1}"
            )
            # add variant edit on the current exonic length
            coding_exonic_length = coding_exonic_length + diff_len
        logger.debug(
            f"after update with variant length, coding_exonic = {coding_exonic_length}, remaining = {remain_codon_length}"
        )
        logger.debug(f"remain codon_length= {remain_codon_length}")
        coding_exonic_length = coding_exonic_length + remain_codon_length
        remain_codon_length = coding_exonic_length % 3
        current_codon_length = coding_exonic_length - remain_codon_length
        logger.debug(
            f"coding exonic length={coding_exonic_length}, remain codon length={remain_codon_length}, current codon length={current_codon_length}"
        )
        # deactivate flag for next exon with coding sequence
        exon_contains_start_codon = False
    elif exon_contains_coding_seq:
        # exon containing coding sequence, but not the start codon
        last_start = last_start + current_codon_length
        logger.debug("updated last start: {last_start}")
        coding_exonic_length = (exon_end - exon_start) + 1
        if coding_exonic_length < 0:
            logger.debug(
                "exon end is lower than exon start, stop codon between two exons"
            )
            coding_pos = ref_transcript.coding_sequence_position_ranges[-1]
            logger.debug(f"coding positions: {coding_pos}")
            coding_exonic_length = (coding_pos[1] - coding_pos[0]) + 1

        logger.debug(
            f"before: adding variant edit, coding_exonic_length={coding_exonic_length}"
        )
        if exon_idx + 1 == exons_containing_var[0] and diff_len != 0:
            logger.debug(
                f"add variant edit len={diff_len} on the current length of exon idx:{exon_idx + 1}"
            )
            coding_exonic_length = coding_exonic_length + diff_len

        logger.debug(
            f"after update with variant length, coding_exonic = {coding_exonic_length}, remaining = {remain_codon_length}"
        )
        coding_exonic_length = coding_exonic_length + remain_codon_length
        remain_codon_length = coding_exonic_length % 3
        current_codon_length = coding_exonic_length - remain_codon_length
        logger.debug(
            f"coding exonic length={coding_exonic_length}, remain codon length={remain_codon_length}, current codon length={current_codon_length}"
        )

    return (
        exon_contains_coding_seq,
        exon_contains_start_codon,
        last_start,
        current_codon_length,
        remain_codon_length,
    )
