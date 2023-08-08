#!/usr/bin/env python3

import logging

import pyensembl
import hgvs.parser
import hgvs.posedit

from refactoring.variant import VariantInfo, TranscriptInfo

logger = logging.getLogger("GenOtoScope_Classify.PVS1.exon_skipping")
hgvs_parser = hgvs.parser.Parser()


def assess_exon_skipping(
    transcript: TranscriptInfo,
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> tuple:
    """
    Assess if exon will be skipped
    By examining if the affected intron position in on the splice sites: +/- 1,2
    """

    logger.debug("Assess exon skipping for splice acceptor variant")

    are_exons_skipped = False
    exons_skipped = []
    start_codon_exon_skipped, stop_codon_exon_skipped = False, False
    coding_exon_skipped = False
    skipped_exon_start, skipped_exon_end = 0, 0

    split_symbols, intron_offsets, directions2exon = parse_variant_intron_pos(
        transcript.var_hgvs
    )

    logger.debug(f"intron offsets: {intron_offsets}")

    for intron_offset in intron_offsets:
        if intron_offset in [1, 2]:
            # variant is disrupting donor/acceptor splicesome site
            # predict that exon will be skipped
            are_exons_skipped = True
    logger.debug(f"is_exon_skipped: {are_exons_skipped}")

    if are_exons_skipped:
        logger.debug("exon is skipped, call find_exon_by_var_pos()")
        # if exon is skipped find its start and end
        exons_skipped, var_exon_start, var_exon_end = find_exon_by_var_pos(
            ref_transcript, transcript, variant, is_genomic=True, diff_len=0
        )
        ### ### ### ###
        # Examine if skipped exon is coding
        ### ### ### ###
        # find exons containing start and stop codons
        if is_transcript_in_positive_strand(ref_transcript):
            start_codon_first_pos = ref_transcript.start_codon_positions[0]
        else:
            start_codon_first_pos = ref_transcript.start_codon_positions[-1]
        start_codon_exon_idx, start_codon_exon_offset = find_exon_by_ref_pos(
            ref_transcript, start_codon_first_pos, True
        )
        stop_codon_exon_idx, stop_codon_exon_offset = find_exon_by_ref_pos(
            ref_transcript, len(ref_transcript.coding_sequence), False
        )

        # The first position of the skipped_exon list is the first exon that is found that overlaps the variant
        if exons_skipped[0] >= start_codon_exon_idx + 1:
            ### ### ### ###
            # skipped exon is coding
            # find start and end position of skipped coding exons
            ### ### ### ###
            coding_exon_skipped = True
            # get as start and end the skipped exon positions
            # AL# List of lists of ints containing the start and end position of the exon
            exon_offsets = get_transcript_exon_offsets(ref_transcript, False)
            try:
                assert len(exons_skipped) == 1
            except AssertionError:
                logger.error(
                    f"Currently intron variant case is implemented only affecting one exon\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )

            # exons offsets contain the length of the sequence of the first exons up to start codon (ATG)
            # so subtract this length to interact with pyensembl transcipt.coding_sequence
            len_first_exons_up_start_codon = ref_transcript.start_codon_spliced_offsets[
                0
            ]
            logger.debug(
                f"len of first exons up to start codon: {len_first_exons_up_start_codon}"
            )
            # AL# Get start and end position of skipped exon
            skipped_exon_offsets = exon_offsets[exons_skipped[0] - 1]
            # AL# Exon start position - distance between exon start and start codon location
            skipped_exon_start = (
                skipped_exon_offsets[0] - len_first_exons_up_start_codon
            )

            if skipped_exon_start <= 1:
                # if exon start is before the start codon position,
                # make the variant start position equal to 1
                skipped_exon_start = 1
                start_codon_exon_skipped = True
            skipped_exon_end = skipped_exon_offsets[1] - len_first_exons_up_start_codon
            logger.debug(
                f"skipped exon start: {skipped_exon_start}, end:{skipped_exon_end}"
            )

            ### ### ###
            # examine if the last exon that was skipped,
            # contained the stop codon
            ### ### ###
            if are_exons_skipped and exons_skipped[0] == stop_codon_exon_idx + 1:
                logger.debug("Search for termination codon on the skipped exon")
                # create skipped coding sequence with in-frame start
                inframe_start = (skipped_exon_start - 1) % 3
                skipped_inframe_seq = ref_transcript.coding_sequence[
                    skipped_exon_start - 1 + inframe_start : skipped_exon_end
                ]
                logger.debug(f"skipped inframe coding seq: {skipped_inframe_seq}")
                if search_termination_codon(extract_codons(skipped_inframe_seq), False):
                    # for stop codon exon skipping, normalize exon end to stop codon position
                    stop_codon_exon_skipped = True
                    skipped_exon_end = len(str(ref_transcript.coding_sequence))
    else:
        logger.debug("Exon is not skipped")
    return (
        exons_skipped,
        are_exons_skipped,
        skipped_exon_start,
        skipped_exon_end,
        start_codon_exon_skipped,
        stop_codon_exon_skipped,
        coding_exon_skipped,
    )


def search_termination_codon(
    exon_codons: list[str], is_penultimate_exon: bool
) -> tuple[bool, int]:
    """
    Search termination codon in exon codons
    Return if termination codon was found and if so where
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


def find_exon_by_ref_pos(
    ref_transcript: pyensembl.transcript.Transcript, ref_pos: int, is_genomic: bool
) -> tuple[int, int]:
    """
    Find exon info that contains the reference position
    the returned exon index is 0-based index
    exon index (0-based) containing the reference position
    offset of reference position in overlapping exon
    """

    logger.debug(f"Find exon containing reference position: {ref_pos}")
    if is_transcript_in_positive_strand(ref_transcript):
        strand_direction = +1
    else:
        strand_direction = -1
    exon_positions = get_transcript_exon_offsets(ref_transcript, is_genomic)
    overlap_exon_index, pos_in_overlap_exon = -1, -1

    # reference position is found using pyensembl coding sequence,
    # this coding sequence does not contain the untranslated region of the starting exons
    # thus you need to add the length of the sequence of the first exons up to start codon (ATG)
    if not is_genomic:
        # for cDNA position don't need to reverse exon start and end
        strand_direction = +1
        len_first_exons_up_start_codon = ref_transcript.start_codon_spliced_offsets[0]
        ref_pos = ref_pos + len_first_exons_up_start_codon

    # find reference position into exon intervals
    for exon_idx, exon_interval in enumerate(exon_positions):
        normalized_exon_interval = range(
            exon_interval[0], exon_interval[1] + strand_direction, strand_direction
        )
        if ref_pos in normalized_exon_interval:
            overlap_exon_index = exon_idx
            pos_in_overlap_exon = abs(ref_pos - exon_interval[0])

    try:
        assert overlap_exon_index >= 0 and pos_in_overlap_exon >= 0
    except AssertionError:
        logger.error(
            f"Exon containing the reference position could not be found\n=> reference pos: {ref_pos}",
            exc_info=True,
        )
    return (overlap_exon_index, pos_in_overlap_exon)


def get_transcript_exon_offsets(
    ref_transcript: pyensembl.transcript.Transcript, is_genomic: bool
) -> list[list[int]]:
    """
    Get all exons cDNA or genomic offsets of transcript
    Return all exons start and end position
    """

    if is_genomic:
        ### ### ###
        # follow strand direction for chromosomal positions:
        # exon_i start < exon_i end for positive strand
        # exon_i start > exon_i end for negative strand
        ### ### ###
        coding_exons_chrom_ranges = []
        for exon in ref_transcript.exons:
            if is_transcript_in_positive_strand(ref_transcript):
                coding_exons_chrom_ranges.append(
                    [exon.to_dict()["start"], exon.to_dict()["end"]]
                )
            else:
                coding_exons_chrom_ranges.append(
                    [exon.to_dict()["end"], exon.to_dict()["start"]]
                )
        return coding_exons_chrom_ranges
    else:
        # create cDNA exon positions
        # offsets starting from 1
        exon_positions = []
        is_first_exon = True
        for exon in ref_transcript.exons:
            if is_first_exon:
                start = 1
                is_first_exon = False
            else:
                start = end + 1
            end = start + exon.to_dict()["end"] - exon.to_dict()["start"]
            exon_positions.append([start, end])
            try:
                assert end - start == exon.to_dict()["end"] - exon.to_dict()["start"]
            except AssertionError:
                logger.error(
                    "Exon genomic regions and coding regions should be equal in size",
                    exc_info=True,
                )
        # sum = 0
        # for exon_position in exon_positions:
        #   sum += exon_position[1] - exon_position[0] + 1
        logger.debug(f"Exon positions: {exon_positions}")
        return exon_positions


def is_transcript_in_positive_strand(
    ref_transcript: pyensembl.transcript.Transcript,
) -> bool:
    """
    Check if transcript is on positive strand
    """
    strand = ref_transcript.exons[0].to_dict()["strand"]
    if strand == "+":
        return True
    else:
        return False


def extract_codons(sequence: str) -> list[str]:
    """
    Extract codons from given sequence, which start where the frame start offset is
    """
    codons = []
    for pos in range(0, len(sequence), 3):
        codons.append(sequence[pos : pos + 3])
    return codons


def is_transcript_type_splice_acceptor_donor(transcript_typ: str) -> bool:
    """
    Examine if transcript type is splice acceptor or splice donor
    """
    intron_variant_types = {
        "splice_acceptor_variant",
        "splice_acceptor",
        "splice_donor_variant",
        "splice_donor",
    }
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

    if (
        transcript_typ in intron_variant_types
        and not transcript_typ in exon_variant_types
    ):
        return True
    else:
        return False


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


def parse_variant_intron_pos(var_coding: hgvs.posedit.PosEdit):
    """
    Parse variant offset in intron

    Returns
    -------
    list of str
        for each variant edit part (start and end):
        '+' to show that variant is just after an exon, '-' to show that the variant is before an exon
    list of int
        for each variant edit part:
        variant offset in intron
    list of int
        for each variant edit part:
        +1 if variant needs to increase its position to get to the closest exon (split_symbol = -)
        -1 if variant needs to decrease its position to get to the closest exon (split_symbol = +)
    """

    logger.debug("Parse intron variant position")
    var_coding_str = str(var_coding)
    var_edit = str(var_coding.edit)
    intron_offset, direction2closest_exon = 0, 0
    split_symbols, intron_offsets, directions2exon = [], [], []
    # find the direction of the closest exon
    if "_" in var_coding_str:
        # for duplication, insertion and deletion
        # split on '_' character before finding direction
        [edit_start, edit_end] = var_coding_str.split(var_edit)[0].split("_")
        for edit_part in var_coding_str.split(var_edit)[0].split("_"):
            if "+" in edit_part:  # after exon in start of the intron
                split_symbol = "+"
                direction2closest_exon = -1
            elif "-" in edit_part:  # before exon in the end of the intron
                split_symbol = "-"
                direction2closest_exon = +1
            else:
                # print("var edit part: {} does not contain split symbol".format(edit_part))
                continue
            # get position of edit part inside the exon

            intron_offset_pos = int(edit_part.split(split_symbol)[1])
            split_symbols.append(split_symbol)
            directions2exon.append(direction2closest_exon)
            intron_offsets.append(intron_offset_pos)
    else:
        # for SNP find direction symbol
        if "+" in var_coding_str:  # after exon in start of the intron
            split_symbol = "+"
            direction2closest_exon = -1
        else:  # before exon in the end of the intron
            split_symbol = "-"
            direction2closest_exon = +1
        split_symbols.append(split_symbol)
        directions2exon.append(direction2closest_exon)
        # to get the splice site position, after splitting on - or +,
        # split on the edit to get the splice offset
        intron_offset_pos = int(
            var_coding_str.split(split_symbol)[-1].split(var_edit)[0]
        )
        intron_offsets.append(intron_offset_pos)
    logger.debug(
        f"split_symbols: {split_symbols}, intron_offsets: {intron_offset}, directions2exon: {directions2exon}".format(
            split_symbols, intron_offsets, directions2exon
        )
    )
    try:
        assert len(split_symbols) == len(intron_offsets) == len(directions2exon)
    except AssertionError:
        logger.error(
            "All parts of variant edit should contain split symbol, intron offset and direction to closest exon.",
            exc_info=True,
        )
    return split_symbols, intron_offsets, directions2exon


def find_exon_by_var_pos(
    ref_transcript: pyensembl.transcript.Transcript,
    transcript: TranscriptInfo,
    variant: VariantInfo,
    is_genomic: bool,
    diff_len: int,
) -> tuple:
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
    # save the variant start and end offset in the overlapping exon regions
    var_start_offset, var_end_offset = (
        0,
        0,
    )
    var_coding = transcript.var_hgvs

    # get the strand direction
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

            if transcript.exon:
                exon_idx = transcript.exon - 1
            elif transcript.intron:
                intron_idx = transcript.intron - 1
                if directions2exon[0] == 1:
                    exon_idx = intron_idx + 1
                else:
                    exon_idx = intron_idx
            else:
                try:
                    assert transcript.exon or transcript.intron
                except AssertionError:
                    logger.error(
                        f"Variant does not contain exon or intron index in VEP column\n => variant position: {variant.to_string()}",
                        exc_info=True,
                    )
            var_start = exon_positions[exon_idx][0]
            var_end = exon_positions[exon_idx][1]
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
            if diff_len > 0:
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
            if var_start in normalized_exon_interval:
                overlap_exon_indices.append(exon_idx + 1)
                break
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
