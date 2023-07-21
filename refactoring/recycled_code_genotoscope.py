#!/usr/bin/env python3

import logging

logger = logging.getLogger("GenOtoScope_Classify.PVS1")

def assess_exon_skipping(transcript, transcript_info, intron_offsets, variant_info):
    """
    Assess if exon will be skipped
    By examining if the affected intron position in on the splice sites: +/- 1,2

    Parameters
    ----------
    transcript : pyensembl.transcript
        ensembl object for variant-affected transcript
    transcript_info : dict of str: str
        variant-affected transcript GSvar information
    intron_offsets : list of int
        variant offset positions in intron
    variant_info : VariantInfo
    variant basic info

    Returns
    -------
    list of int
        list of exon (1-based) indices containing the variant
    bool
        True if exon predicted to be skipped, otherwise False
    int
        if exon is skipped, variant deletion start position
    int
        if exon is skipped, variant deletion end position
    bool
        variant skips start codon exon
    bool
        variant skips stop codon exon
    bool
        variant skips coding exon
    """

    print("Assess exon skipping for splice acceptor variant")
    is_exon_skipped = False
    variant_skips_start_codon_exon, variant_skips_stop_codon_exon = False, False
    variant_skips_coding_exon = False
    skipped_exon_start, byiskipped_exon_end = 0, 0
    print("intron offsets: {}".format(intron_offsets))

    for intron_offset in intron_offsets:
        if intron_offset in [1, 2]:
            # variant is disrupting donor/acceptor splicesome site
            # predict that exon will be skipped
            is_exon_skipped = True
    print("is_exon_skipped: {}".format(is_exon_skipped))

    if is_exon_skipped:
        print("exon is skipped, call find_exon_by_var_pos()")
        # if exon is skipped find its start and end
        skipped_exons, var_exon_start, var_exon_end = find_exon_by_var_pos(transcript,
                                                                                transcript_info,
                                                                                variant_info,
                                                                                is_genomic=True)
        ### ### ### ###
        # Examine if skipped exon is coding
        ### ### ### ###
        # find exons containing start and stop codons
        if is_transcript_in_positive_strand(transcript):
            transcript_strand = "+"
            start_codon_first_pos = transcript.start_codon_positions[0]
        else:
            transcript_strand = "-"
            start_codon_first_pos = transcript.start_codon_positions[-1]
        start_codon_exon_idx, start_codon_exon_offset = find_exon_by_ref_pos(transcript, start_codon_first_pos,
                                                                                  True)
        stop_codon_exon_idx, stop_codon_exon_offset = find_exon_by_ref_pos(transcript,
                                                                                len(transcript.coding_sequence),
                                                                                False)
        # The first position of the skipped_exon list is the first exon that is found that overlaps the variant
        if skipped_exons[0] >= start_codon_exon_idx + 1:
            ### ### ### ###
            # skipped exon is coding
            # find start and end position of skipped coding exons
            ### ### ### ###
            variant_skips_coding_exon = True
            # get as start and end the skipped exon positions
            #AL# List of lists of ints containing the start and end position of the exon
            exon_offsets = get_transcript_exon_offsets(transcript, False)
            try:
                assert len(skipped_exons) == 1
            except AssertionError:
               logger.error(
                    "Currently intron variant case is implemented only affecting one exon\n=> variant position: {}".format(
                        variant_info.to_string()), exc_info=True)

            # exons offsets contain the length of the sequence of the first exons up to start codon (ATG)
            # so subtract this length to interact with pyensembl transcipt.coding_sequence
            len_first_exons_up_start_codon = transcript.start_codon_spliced_offsets[0]
            logger.debug("len of first exons up to start codon: {}".format(len_first_exons_up_start_codon))
            #AL# Get start and end position of skipped exon
            skipped_exon_offsets = exon_offsets[skipped_exons[0] - 1]
            #AL# Exon start position - distance between exon start and start codon location
            skipped_exon_start = skipped_exon_offsets[0] - len_first_exons_up_start_codon

            if skipped_exon_start <= 1:
                # if exon start is before the start codon position,
                # make the variant start position equal to 1
                skipped_exon_start = 1
                variant_skips_start_codon_exon = True
            skipped_exon_end = skipped_exon_offsets[1] - len_first_exons_up_start_codon
            logger.debug("skipped exon start: {}, end:{}".format(skipped_exon_start, skipped_exon_end))

            ### ### ###
            # examine if the last exon that was skipped,
            # contained the stop codon
            ### ### ###
            if is_exon_skipped and skipped_exons[0] == stop_codon_exon_idx + 1:
                logger.debug("Search for termination codon on the skipped exon")
                # create skipped coding sequence with in-frame start
                inframe_start = (skipped_exon_start - 1) % 3
                skipped_inframe_seq = transcript.coding_sequence[
                    skipped_exon_start - 1 + inframe_start:skipped_exon_end]
                logger.debug("skipped inframe coding seq: {}".format(skipped_inframe_seq))
                if search_termination_codon(extract_codons(skipped_inframe_seq), False):
                    # for stop codon exon skipping, normalize exon end to stop codon position
                    variant_skips_stop_codon_exon = True
                    skipped_exon_end = len(str(transcript.coding_sequence))
    else:
        logger.debug("Exon is not skipped")
    return skipped_exons, is_exon_skipped, skipped_exon_start, skipped_exon_end, variant_skips_start_codon_exon, variant_skips_stop_codon_exon, variant_skips_coding_exon


def find_exon_by_var_pos(transcript, transcript_info, variant_info, is_genomic):
    """
    Find variant exon index by variant coding position, exon index is 1-based

    Parameters
    ----------
    transcript : pyensembl.transcript
        transcript object contains the penultimate and final exon
    transcript_info : dict of str: str
        transcript information
    variant_info : VariantInfo
        variant basic info
    is_genomic : bool
        genomic positions will be used (True), otherwise cDNA positions will be used (False)

    Returns
    -------
    list of int
        indices of exons overlapping with variant (1-based)
    int
        variant start position as offset in first overlapping exon
    int
        variant end	position as offset in last overlapping exon
    """

    logger.debug("Find exon indices containing the variant")
    if not ("exon" in transcript_info["exon"] or "intron" in transcript_info["exon"]):
        logger.debug("Exon information is not found in VEP => use VEP spliced offset")
        is_genomic = False

    exon_positions = get_transcript_exon_offsets(transcript, is_genomic)

    overlap_exon_indices = []
    var_start_offset, var_end_offset = 0, 0  # save the variant start and end offset in the overlapping exon regions
    var_coding = transcript_info["var_coding"]
    ### ### ###
    # get the strand direction
    ### ### ###
    if is_transcript_in_positive_strand(transcript):
        strand_direction = + 1
    else:
        strand_direction = - 1
    if is_genomic:
        if not is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
            ### ### ###
            # Exonic variant
            ### ### ###
            logger.debug("Exonic variant type")
            var_start = variant_info.genomic_start
            var_end = variant_info.genomic_end
        else:
            logger.debug("Intron variant type => update as variant pos the starting position of skipped exon")
            split_symbols, intron_offsets, directions2exon = parse_variant_intron_pos(var_coding)

            ### ### ###
            # Intronic variant
            # Parse affected exon index by VEP (1-based)
            # add splice site direction to affected exon index
            # (currently modeled that) intron variant can disrupt at most one exon, so equal variant end with its start position
            ### ### ###

            if "exon" in transcript_info["exon"] and "intron" in transcript_info["exon"]:
                exon_idx = int(transcript_info["exon"].split("intron")[1].split("/")[0]) - 1
                logger.debug("Transcript info contains both exon and intron offset => use intron offset")
            elif "exon" in transcript_info["exon"]:
                exon_idx = int(transcript_info["exon"].split("/")[0].split("exon")[1]) - 1
            elif "intron" in transcript_info["exon"]:
                exon_idx = int(transcript_info["exon"].split("/")[0].split("intron")[1]) - 1
            else:
                try:
                    assert "exon" in transcript_info["exon"] or "intron" in transcript_info["exon"]
                except AssertionError:
                    logger.error(
                        "Variant does not contain exon or intron index in VEP column\n => variant position: {}".format(
                            variant_info.to_string()), exc_info=True)
            if exon_idx + directions2exon[0] >= 0:
                skipped_exon = exon_positions[exon_idx + directions2exon[0]]
            else:
                logger.debug("Skipping first exon")
                skipped_exon = exon_positions[0]
                var_start = skipped_exon[0]
                var_end = var_start
                logger.debug("Updated variant start:{}, end: {} on exon idx: {}".format(var_start, var_end,
                                                                                             exon_idx + directions2exon[
                                                                                                 0]))
    else:
        # use cDNA offset for both frameshift and nonsense mutation
        var_start = int(var_coding.pos.start.base)
        var_end = int(var_coding.pos.end.base)
        if var_start < 0:
            var_start = 0
        if var_end < 0:
            if transcript_info["coding_diff_len"] > 0:
                var_end = var_start + transcript_info["coding_diff_len"]
            else:
                # deletion case
                var_end = var_start

        # VEP positions include only the coding sequence of the transcript,
        # so you need to add the length of the sequence of the first exons up to start codon (ATG)
        len_first_exons_up_start_codon = transcript.start_codon_spliced_offsets[0]
        var_start = var_start + len_first_exons_up_start_codon
        var_end = var_end + len_first_exons_up_start_codon

    ### ### ###
    # find variant position into exon intervals
    ### ### ###
    for exon_idx, exon_interval in enumerate(exon_positions):
        if is_genomic:
            normalized_exon_interval = range(exon_interval[0], exon_interval[1] + 2 * strand_direction,
                                             strand_direction)
        else:
            normalized_exon_interval = range(exon_interval[0], exon_interval[1] + 1)
            logger.debug("Exon interval: {}".format(normalized_exon_interval))
            if var_start in normalized_exon_interval:
                overlap_exon_indices.append(exon_idx + 1)
                break
    logger.debug("Search var_start: {} found in exon(s): {}".format(var_start, overlap_exon_indices))

    ### ### ###
    # find variant start and end offset in exons
    ### ### ###
    if len(overlap_exon_indices) == 1:  # variant included in only one exon
        if is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
            # exon-skipping, get as offset the total range of skipped exon
            var_start_offset = 0
            # negative strand transcripts contain higher start position than end position
            var_end_offset = abs(exon_positions[overlap_exon_indices[0] - 1][1] - \
                                 exon_positions[overlap_exon_indices[0] - 1][0])
        else:
            if is_genomic:
                if strand_direction == 1:
                    # for positive strand, compute distance from exon lower position (start)
                    var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
                    var_end_offset = var_end - exon_positions[overlap_exon_indices[0] - 1][0]
                else:
                    # for negative strand, compute distance from exon higher position (end)
                    var_start_offset = exon_positions[overlap_exon_indices[0] - 1][0] - var_end
                    var_end_offset = exon_positions[overlap_exon_indices[0] - 1][0] - var_start
            else:
                var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
                var_end_offset = var_end - exon_positions[overlap_exon_indices[0] - 1][0]
    else:
        try:
            assert len(overlap_exon_indices) > 0
        except AssertionError:
            logger.error("Overlapping exon indices should be more than 0\n=> variant position: {}".format(
                variant_info.to_string()), exc_info=True)
        ### ### ###
        # variant included in more than one exons,
        # so start is in the first exon, end is in the last exon
        ### ### ###
        if is_genomic:
            if strand_direction == 1:
                var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
                var_end_offset = var_end - exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][0]
            else:
                # for negative strand, use the highest value (exon start) to compute the offset for the start and exon position of the variant
                var_start_offset = exon_positions[overlap_exon_indices[0] - 1][0] - var_end
                var_end_offset = exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][0] - var_start
        else:
            var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
            var_end_offset = var_end - exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][0]
    try:
        assert var_start_offset >= 0 and var_end_offset >= 0
    except AssertionError:
        logger.error("Variant start and end offset should be higher than 0\n=> variant position: {}".format(
            variant_info.to_string()), exc_info=True)
        logger.debug(
            "Overlap exon indices: {}, var start offset: {}, var end offset: {}".format(overlap_exon_indices,
                                                                                        var_start_offset,
                                                                                        var_end_offset))
    return overlap_exon_indices, var_start_offset, var_end_offset
