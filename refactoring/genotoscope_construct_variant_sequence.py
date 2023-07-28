#!/usr/bin/env python3

import logging
import pyensembl

from genotoscope_exon_skipping import (
    is_transcript_type_splice_acceptor_donor,
)
from genotoscope_assess_NMD import is_genomic_pos_in_coding_exon
from variant import TranscriptInfo, VariantInfo

logger = logging.getLogger("GenOtoScope_Classify.PVS1.construct_variant_coding_seq")


def construct_variant_coding_seq_intronic_variant(
    transcript: TranscriptInfo,
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    var_start: int,
    var_end: int,
    are_exons_skipped: bool,
    start_codon_exon_skipped: bool,
    stop_codon_exon_skipped: bool,
    coding_exon_skipped: bool,
) -> tuple:
    """
    Construct coding sequence for intronic variant
    following hgvs recommendations: http://varnomen.hgvs.org/recommendations/general/

    Then return variant-integrated coding sequence in the forward strand
    r the reverse complement, based on the transcript directionality

    Returns
    -------
    str
        variant coding sequence
    int
        diff_len
    """
    logger.debug("Constructing coding sequence for intronic variant")
    ref_coding_seq = str(ref_transcript.coding_sequence)
    logger.debug(f"Coding seq: {ref_coding_seq}, len={len(ref_coding_seq)}")
    var_coding_seq = ""
    diff_len = 0
    if start_codon_exon_skipped and stop_codon_exon_skipped:
        ### ### ### ###
        # variant skips both exon containing both start and stop codon
        ### ### ### ###
        var_coding_seq = ""
        diff_len = len(ref_coding_seq)
    elif coding_exon_skipped:
        ### ### ### ###
        # Simulate coding exon skipping
        # by deletion of coding exon sequence
        ### ### ### ###
        if are_exons_skipped:
            # exon is skipped => delete exon from coding sequence
            logger.debug("Delete the skipped exon from the coding sequence")
            if var_start >= 1:
                var_coding_seq = (
                    ref_coding_seq[0 : var_start - 1]
                    + ref_coding_seq[var_end : len(ref_coding_seq)]
                )
                logger.debug(
                    f"Deleting sequence: {ref_coding_seq[var_start -1 : var_end]}"
                )
            else:
                try:
                    assert var_start >= 1
                except AssertionError:
                    logger.error(
                        f"Supplied coding region position should be >= 0\n=> variant position: {variant.to_string()}",
                        exc_info=True,
                    )

            # calculate deletion length and assert its result
            del_length = var_end - var_start + 1
            logger.debug(
                f"len(coding)={len(ref_coding_seq)}, len(var_coding)={len(var_coding_seq)}, del_length={del_length}"
            )
            try:
                assert len(ref_coding_seq) - len(var_coding_seq) == del_length
            except AssertionError:
                logger.error(
                    f"For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n=> variant position: {variant.to_string}",
                    exc_info=True,
                )
        diff_len = -1 * del_length
    else:
        var_coding_seq = ref_coding_seq

    if coding_exon_skipped:
        try:
            assert ref_coding_seq.upper() != var_coding_seq.upper()
        except AssertionError:
            logger.error(
                f"Coding sequence of reference and sample should be different\n variant position: {variant.to_string()}",
                exc_info=True,
            )
        print_ref_observed_seq(
            ref_coding_seq,
            var_coding_seq,
            transcript,
            var_start,
            var_end,
            var_edit,
            are_exons_skipped,
        )

        ### ### ### ### ### ###
        # assert constructed variant coding sequence
        # contains start (if edit not start_lost) and stop codon
        ### ### ### ### ### ###
        if start_codon_exon_skipped and stop_codon_exon_skipped:
            pass
        elif start_codon_exon_skipped:
            try:
                assert var_coding_seq[-3:] in ["TAG", "TAA", "TGA"]
            except AssertionError:
                logger.error(
                    f"Constructed observed coding sequence should contain stop codon\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )
        elif stop_codon_exon_skipped:
            if ref_coding_seq[0:3] == "ATG":
                try:
                    assert var_coding_seq[0:3] == "ATG"
                except AssertionError:
                    logger.debug(
                        f"Constructed observed coding sequence should contain start codon\n=> variant position: {variant.to_string()}"
                    )
        else:
            check_if_start_stop_codon_exists_in_var_sequence(
                transcript, variant, var_start, var_end, ref_coding_seq, var_coding_seq
            )
    return var_coding_seq, diff_len


def construct_variant_coding_seq_exonic_variant(
    transcript: TranscriptInfo,
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> tuple:
    """
    Construct coding sequence for exonic variant
    following hgvs recommendations: http://varnomen.hgvs.org/recommendations/general/

    Then return variant-integrated coding sequence in the forward strand
    r the reverse complement, based on the transcript directionality

    Returns
    -------
    str
        variant coding sequence
    int
        diff_len
    """
    logger.debug("Constructing coding sequence for intronic variant")
    ref_coding_seq = str(ref_transcript.coding_sequence)
    logger.debug(f"Coding seq: {ref_coding_seq}, len={len(ref_coding_seq)}")
    vep_coordinates_used = True
    var_coding_seq = ""
    diff_len = 0
    var_edit = str(transcript.var_hgvs.edit)
    logger.debug("Constructing variant start and stop position")
    if int(transcript.var_hgvs.pos.start.base) >= 1:
        var_start = int(transcript.var_hgvs.pos.start.base)
    else:
        var_start = 1
    if int(transcript.var_hgvs.pos.end.base) >= 1:
        var_end = int(transcript.var_hgvs.pos.end.base)
    else:
        var_end = 1

    ### ### ### ###
    # construct observed coding sequence for variant in exonic regions
    ### ### ### ###
    if ">" in var_edit[0:3]:
        logger.debug("Add SNP to coding sequence")
        # get directionally corrected reference and observed bases for SNP
        [ref_seq, obs_seq] = var_edit.split(">")
        if var_start > 1:
            var_coding_seq = (
                ref_coding_seq[0 : var_start - 1]
                + obs_seq.lower()
                + ref_coding_seq[var_start : len(ref_coding_seq)]
            )
        elif var_start == 1:
            var_coding_seq = (
                obs_seq.lower() + ref_coding_seq[var_start : len(ref_coding_seq)]
            )
        else:
            try:
                assert var_start >= 0
            except AssertionError:
                logger.error(
                    f"Supplied coding region position should be >= 0\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )
        try:
            assert len(ref_coding_seq) == len(var_coding_seq)
        except AssertionError:
            logger.error(
                f"For SNP the sum of length of coding exons should not change\n=> variant position: {variant.to_string()}",
                exc_info=True,
            )
        diff_len = 0
    elif "delins" in var_edit[0:6]:
        ### ### ###
        # delins variant edits add a deletion and then an insertion
        ### ### ###
        logger.debug("Add deletion and then insertion to coding sequence")

        ### ### ###
        # introduce the deletion
        ### ### ###
        logger.debug(f"Deleting: {ref_coding_seq[var_start -1 :var_end]}")
        logger.debug(f"Region being delete is: {var_start} - {var_end}")
        if var_start >= 1:
            var_coding_seq = (
                ref_coding_seq[0 : var_start - 1]
                + ref_coding_seq[var_end : len(ref_coding_seq)]
            )
        elif var_start == 0:
            var_coding_seq = ref_coding_seq[var_end : len(ref_coding_seq)]
        else:
            try:
                assert var_start >= 0
            except AssertionError:
                logger.error(
                    f"Supplied coding region position should be >= 0\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )

        # assert deletion by resulted sequence length
        del_length = var_end - var_start + 1
        try:
            logger.debug(
                f"var_coding_seq: {var_coding_seq}, len: {len(var_coding_seq)}"
            )
            assert len(ref_coding_seq) - len(var_coding_seq) == del_length
        except AssertionError:
            logger.error(
                f"For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n variant position: {variant.to_string()}",
                exc_info=True,
            )
        ### ### ###
        # introduce the insertion
        ### ### ###
        obs_seq = var_edit.split("delins")[1]
        try:
            assert len(obs_seq) >= 1
        except AssertionError:
            logger.error(
                f"Supplied VEP does not contain insertion sequence\n=> variant position: {variant.to_string()}",
                exc_info=True,
            )
        logger.debug(f"Insert sequence: {obs_seq}")
        if var_start >= 1:
            var_coding_seq_ins = (
                var_coding_seq[0 : var_start - 1]
                + obs_seq.lower()
                + var_coding_seq[var_start - 1 : len(ref_coding_seq)]
            )
        else:
            try:
                assert var_start >= 1
            except AssertionError:
                logger.error(
                    f"Supplied coding region position should be >= 1\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )

        # assert insertion operation by resulted length
        ins_length = len(obs_seq)
        try:
            assert len(var_coding_seq_ins) - len(var_coding_seq) == ins_length
        except AssertionError:
            logger.error(
                f"For insertion should hold: len(var_coding) - len(var_coding with deletion) = len(insertion)\n=> variant position: {variant.to_string()}",
                exc_info=True,
            )
        var_coding_seq = var_coding_seq_ins
        # assert the length difference for both operations
        diff_len = ins_length - del_length
        try:
            logger.debug(
                f"var_coding_seq: {var_coding_seq}, len: {len(var_coding_seq)}"
            )
            assert len(var_coding_seq) - len(ref_coding_seq) == diff_len
        except AssertionError:
            logger.error(
                f"For insertion after deletion should hold: len(variant coding) - len(reference coding) = -len(deletion) + len(insertion)\n variant position: {variant.to_string()}",
                exc_info=True,
            )
    elif "del" in var_edit[0:3]:
        # deletion will be from the start up to end position
        logger.debug("Add deletion to coding sequence")
        if not is_transcript_type_splice_acceptor_donor(transcript.var_type):
            # assert deleted coding sequence to equal the reference coding sequence
            # include case at which the variant start before or ends after the coding sequence
            if (
                int(transcript.var_hgvs.pos.start.base) > 0
                and int(transcript.var_hgvs.pos.end.base) > 0
            ):
                if not (
                    is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_start)
                    and is_genomic_pos_in_coding_exon(
                        ref_transcript, variant.genomic_end
                    )
                ):
                    # variant start or ends outside an coding exonic range
                    logger.debug("Variant start or end not in a coding exonic range")
                    logger.debug(f"var start: {var_start}, end: {var_end}")
                    affected_coding_length = var_end - var_start + 1
                    logger.debug(
                        f"Length of variant sequence that affects coding exonic sequence: {affected_coding_length}"
                    )
                else:
                    # variant starts and ends inside coding exonic range
                    affected_coding_length = -1
        ### ### ###
        # perform deletion
        ### ### ###
        if var_start >= 1:
            var_coding_seq = (
                ref_coding_seq[0 : var_start - 1]
                + ref_coding_seq[var_end : len(ref_coding_seq)]
            )
        else:
            try:
                assert var_start >= 0
            except AssertionError:
                logger.error(
                    f"Supplied coding region position should be >= 0\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )
        del_length = var_end - var_start + 1
        try:
            logger.debug(
                f"var_coding_seq: {var_coding_seq}, len: {len(var_coding_seq)}"
            )
            assert len(ref_coding_seq) - len(var_coding_seq) == del_length
        except AssertionError:
            logger.error(
                f"For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n variant position: {variant.to_string()}",
                exc_info=True,
            )
        diff_len = -1 * del_length
    elif "ins" in var_edit[0:3]:
        # for insertion annotation: start & end are the flanking regions
        # the insertion will be placed between the flanking regions
        logger.debug("Add insertion to coding sequence")
        if not (
            is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_start)
            and is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_end)
        ):
            # variant start or ends outside an coding exonic range
            logger.debug("Variant start or end not in a coding exonic range")
            affected_coding_length = var_end - var_start + 1
            logger.debug(
                f"Length of variant sequence that affects coding exonic sequence: {affected_coding_length}"
            )
        else:
            # variant inside coding exonic range
            affected_coding_length = -1

        ### ### ###
        # calculate sequence to be inserted
        ### ### ###
        obs_seq = var_edit.split("ins")[1]
        if len(obs_seq) > 0:
            # VEP contains insertion sequence
            if "_" in obs_seq:
                # insertion sequence is described by coding coordinates
                logger.info(
                    f"Insertion sequence is described by coding coordinates \n=> variant pos: {variant.to_string}"
                )
                [ins_source_start, ins_source_end] = obs_seq.split("_")
                obs_seq = ref_coding_seq[
                    int(ins_source_start) - 1 : int(ins_source_end.strip())
                ]
            if not affected_coding_length == -1:
                obs_seq = obs_seq[-affected_coding_length:]
        else:
            try:
                assert 1 == 0
            except AssertionError:
                logger.error(
                    f"vep does not contain insertion sequence\n=> variant position: {variant.to_string}",
                    exc_info=True,
                )
        logger.debug(f"Insert the sequence: {obs_seq}")
        ### ### ###
        # perform insertion
        ### ### ###
        if var_start >= 1:
            if var_start == var_end:
                # start equal ends so the rest part, after the insertion, should start on end position
                if vep_coordinates_used:
                    var_coding_seq = (
                        ref_coding_seq[0:var_start]
                        + obs_seq.lower()
                        + ref_coding_seq[var_end : len(ref_coding_seq)]
                    )
                else:
                    var_coding_seq = (
                        ref_coding_seq[0:var_start]
                        + obs_seq.lower()
                        + ref_coding_seq[var_end : len(ref_coding_seq)]
                    )
            else:
                if vep_coordinates_used:
                    var_coding_seq = (
                        ref_coding_seq[0:var_start]
                        + obs_seq.lower()
                        + ref_coding_seq[var_end - 1 : len(ref_coding_seq)]
                    )
                else:
                    var_coding_seq = (
                        ref_coding_seq[0:var_start]
                        + obs_seq.lower()
                        + ref_coding_seq[var_end - 1 : len(ref_coding_seq)]
                    )
        else:
            try:
                assert var_start >= 1
            except AssertionError:
                logger.error(
                    f"Supplied coding region position should be >= 0\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )
        ins_length = len(obs_seq)
        try:
            assert len(var_coding_seq) - len(ref_coding_seq) == ins_length
        except AssertionError:
            logger.error(
                f"For insertion should hold: len(var_coding) - len(reference_coding) = len(insertion)\n=> variant position: {variant.to_string()}",
                exc_info=True,
            )
        diff_len = +1 * ins_length

    elif "dup" in var_edit[0:3]:
        ### ### ###
        # for duplication annotation: start & end are the duplication region
        # the duplication will be placed right-after the end position
        ### ### ###
        logger.debug("Add duplication to coding sequence")
        if not (
            is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_start)
            and is_genomic_pos_in_coding_exon(ref_transcript, variant.genomic_end)
        ):
            ### ### ###
            # variant starts or ends in intronic region
            ### ### ###
            logger.debug("Variant start or end not in a coding exonic range")
            logger.debug(f"var start: {var_start}, end: {var_end}")
            affected_coding_length = var_end - var_start + 1
            logger.debug(
                f"Length of variant sequence that affects coding exonic sequence: {affected_coding_length}"
            )
        else:
            ### ### ###
            # variant inside coding exonic range
            ### ### ###
            affected_coding_length = -1

        ### ### ###
        # calculate sequence to be duplicated
        ### ### ###
        obs_seq = (
            ("c." + str(transcript.var_hgvs))
            .split("c.")[1]
            .split(str(transcript.var_hgvs))[1]
        )
        if len(obs_seq) > 0:
            if affected_coding_length == -1:
                # all duplication inside coding exon
                obs_seq = obs_seq
            else:
                obs_seq = obs_seq[-affected_coding_length:]
        else:
            if var_start == var_end:
                obs_seq = ref_coding_seq[var_start - 1]
            else:
                obs_seq = ref_coding_seq[var_start - 1 : var_end]
        logger.debug(f"Duplicate the sequence: {obs_seq}")

        ### ### ###
        # perform duplication
        ### ### ###
        if var_start == var_end:
            # duplication of one nucleotide
            logger.debug("Duplication of one nucleotide")
            if var_start >= 1:
                var_coding_seq = (
                    ref_coding_seq[0:var_start]
                    + obs_seq.lower()
                    + ref_coding_seq[var_start : len(ref_coding_seq)]
                )
            else:
                try:
                    assert var_start >= 0
                except AssertionError:
                    logger.error(
                        f"Supplied coding region position should be >= 0\n=> variant position: {variant.to_string()}",
                        exc_info=True,
                    )
        else:
            # duplication of multiple nucleotides
            logger.debug("Duplication of multiple nucleotides")
            if var_start >= 1:
                if vep_coordinates_used:
                    var_coding_seq = (
                        ref_coding_seq[0:var_end]
                        + obs_seq.lower()
                        + ref_coding_seq[var_end : len(ref_coding_seq)]
                    )
                else:
                    var_coding_seq = (
                        ref_coding_seq[0 : var_end + 1]
                        + obs_seq.lower()
                        + ref_coding_seq[var_end + 1 : len(ref_coding_seq)]
                    )
            else:
                try:
                    assert var_start >= 1
                except AssertionError:
                    logger.error(
                        f"Supplied coding region position should be >= 0\n=> variant position: {variant.to_string()}",
                        exc_info=True,
                    )
        dupl_length = len(obs_seq)
        try:
            assert len(var_coding_seq) - len(ref_coding_seq) == dupl_length
        except AssertionError:
            logger.error(
                f"For duplication should hold: len(var_coding) - len(reference_coding) = len(duplication)\n variant position: {variant.to_string()}",
                exc_info=True,
            )
            diff_len = +1 * dupl_length

    ### ### ### ###
    # after creating coding sequence with variant,
    # check that indeed is different from the reference
    ### ### ### ###
    try:
        assert ref_coding_seq.upper() != var_coding_seq.upper()
    except AssertionError:
        logger.error(
            f"Coding sequence of reference and sample should be different\n variant position: {variant.to_string()}",
            exc_info=True,
        )
    print_ref_observed_seq(
        ref_coding_seq,
        var_coding_seq,
        transcript,
        var_start,
        var_end,
        var_edit,
        False,
    )

    ### ### ### ### ### ###
    # assert constructed variant coding sequence
    # contains start (if edit not start_lost) and stop codon
    ### ### ### ### ### ###
    check_if_start_stop_codon_exists_in_var_sequence(
        transcript, variant, var_start, var_end, ref_coding_seq, var_coding_seq
    )
    return var_coding_seq, diff_len


def check_if_start_stop_codon_exists_in_var_sequence(
    transcript: TranscriptInfo,
    variant: VariantInfo,
    var_start: int,
    var_end: int,
    ref_coding_seq: str,
    var_coding_seq: str,
) -> None:
    if var_start == var_end:
        var_positions = range(var_start, var_start + 1)
    else:
        var_positions = range(var_start, var_end)

    if not (
        "start_lost" in transcript.var_type or "start_retained" in variant.var_type
    ):
        start_positions = set(range(1, 4))
        if (
            not start_positions.intersection(var_positions)
            and ref_coding_seq[0:3] == "ATG"
        ):
            # if variant does not change positions on the start of the coding sequence
            # and ensembl record for transcript starts by ATG
            # assert that the observed coding sequence starts with ATG
            try:
                assert var_coding_seq[0:3] == "ATG"
            except AssertionError:
                logger.error(
                    f"Constructed observed coding sequence should contain start codon\n=> variant position: {variant.to_string()}",
                    exc_info=True,
                )

        if not (
            "stop_retained_variant" in transcript.var_type
            or "stop_lost" in transcript.var_type
        ):
            stop_positions = set(
                range(len(ref_coding_seq) + 1 - 3, len(ref_coding_seq) + 1)
            )
            if not stop_positions.intersection(var_positions):
                # if variant does not change positions in the end of the coding sequence
                # assert that the observed coding sequence finishes with "TAG", "TAA" or "TGA"
                try:
                    assert var_coding_seq[-3:] in ["TAG", "TAA", "TGA"]
                except AssertionError:
                    logger.error(
                        f"Constructed observed coding sequence should contain stop codon\n=> variant position: {variant.to_string()}",
                        exc_info=True,
                    )


def print_ref_observed_seq(
    ref_coding_seq: str,
    var_coding_seq: str,
    transcript: TranscriptInfo,
    var_start: int,
    var_end: int,
    var_edit: str,
    is_exon_skipped: bool,
) -> None:
    """
    Print reference and observed coding sequence
    """
    logger.debug("Print reference and observed sequence on the variant region:")
    # set up the variant print region start and end
    if var_start >= 11:
        # for variants after the 11th genomic position, print 10 before and 10 after bases
        print_start = 25
        print_end = 25
    elif 1 < var_start < 11:
        # for variants close to start of the chromosome, print bases as much as to show variant
        # and 10 bases after
        print_start = var_start - 1
        print_end = 25
    else:
        # for variants on the first base of the chromosome, print the start of the chromosome and 10 bases after
        print_start = 0
        print_end = 25
    ### ### ###
    # print variant coding sequence and coding sequence
    ### ### ###
    if "del" in var_edit[0:3] or (
        is_transcript_type_splice_acceptor_donor(transcript.var_type)
        and is_exon_skipped
    ):
        # deletion or skipped exon case
        ref_deleted_var = (
            ref_coding_seq[var_start - print_start - 1 : var_start - 1]
            + ref_coding_seq[var_start - 1 : var_end].lower()
            + ref_coding_seq[var_end : var_end + print_end]
        )
        logger.debug(f">coding seq:\n {ref_deleted_var}")
    else:
        # insertion, duplication, SNP
        logger.debug(
            f">coding seq:\n {ref_coding_seq[var_start - print_start -1 : var_end + print_end]}"
        )

    if var_start == 0:
        logger.debug(
            f">var coding seq:\n {var_coding_seq[var_start : var_end + print_end]}"
        )
    else:
        logger.debug(
            f">var coding seq:\n {var_coding_seq[var_start - print_start -1 : var_end + print_end]}"
        )
