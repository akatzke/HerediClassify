#!/usr/bin/env python3

import logging
import pyensembl

from genotoscope_exon_skipping import (
    get_transcript_exon_offsets,
    is_transcript_in_positive_strand,
    assess_exon_skipping,
    find_exon_by_ref_pos,
    is_transcript_type_splice_acceptor_donor,
    parse_variant_intron_pos,
)
from genotoscope_assess_NMD import is_genomic_pos_in_coding_exon
from genotoscope_refactored_unused import find_exon_by_var_pos
from refactoring.variant import TranscriptInfo, VariantInfo

logger = logging.getLogger("GenOtoScope_Classify.PVS1.construct_variant_coding_seq")


def construct_variant_coding_seq(transcript, transcript_info, variant_info):
    """
    Add variant to coding sequence to get observed coding sequence
    following hgvs recommendations: http://varnomen.hgvs.org/recommendations/general/

    Then return variant-integrated coding sequence in the forward strand
    r the reverse complement, based on the transcript directionality

    Parameters
    ----------
    transcript : pyensembl.transcript
        ensembl object for variant-affected transcript
    transcript_info : dict of str: str
        variant-affected transcript GSvar information
    variant_info : VariantInfo
        variant basic info

    Returns
    -------
    str
        variant coding sequence
    int
        diff_len
    """

    logger.debug(
        "Add variant to coding sequence to create observed coding sequence of transcript id: {}".format(
            transcript_info["transcript_id"]
        )
    )
    var_coding = transcript_info["var_coding"]
    var_edit = str(var_coding.edit)

    is_exon_skipped = False
    diff_len = 0
    (
        variant_skips_start_codon_exon,
        variant_skips_stop_codon_exon,
        variant_skips_coding_exon,
    ) = (False, False, False)
    # find start end of variant positions
    if not is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
        logger.debug("Exonic variant")
        # subtract from variant cDNA position the length of the prime 5' end,
        # to get position in only coding parts of exons
        if int(var_coding.pos.start.base) >= 1:
            var_start = int(var_coding.pos.start.base)
        else:
            var_start = 1
        if int(var_coding.pos.end.base) >= 1:
            var_end = int(var_coding.pos.end.base)
        else:
            var_end = 1
    else:
        logger.debug("Intronic variant")
        # find if variant results to exon skipping
        (
            exons_containing_var,
            is_exon_skipped,
            var_start,
            var_end,
            variant_skips_start_codon_exon,
            variant_skips_stop_codon_exon,
            variant_skips_coding_exon,
        ) = assess_exon_skipping(transcript, transcript_info, variant_info)
        try:
            assert is_exon_skipped is True
        except AssertionError:
            logger.error(
                "Creating variant observed coding sequence, implemented for exon skipping case only\n=> variant position: {}".format(
                    variant_info.to_string()
                ),
                exc_info=True,
            )
        logger.debug(
            "Exon(s) in list: {} are skipped".format(
                exons_containing_var, is_exon_skipped
            )
        )
        logger.debug(
            "variant skips start codon: {}, stop codon: {}, skips coding exon: {}".format(
                variant_skips_start_codon_exon,
                variant_skips_stop_codon_exon,
                variant_skips_coding_exon,
            )
        )
        logger.debug(
            "Skipped exon offset start: {}, end: {}".format(var_start, var_end)
        )
    # load the coding sequence to introduce the variant
    coding_seq = str(transcript.coding_sequence)
    logger.debug("Coding seq: {}, len={}".format(coding_seq, len(coding_seq)))
    var_coding_seq = ""
    vep_coordinates_used, vep_contains_seq = True, False

    ### ### ### ###
    # integrate variant into coding sequence, per variant type
    ### ### ### ###
    if is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
        # to delete skipped exon sequence
        # use the calculated ensembl coding offsets
        vep_coordinates_used, vep_contains_seq = False, False
        if variant_skips_start_codon_exon and variant_skips_stop_codon_exon:
            ### ### ### ###
            # variant skips both exon containing both start and stop codon
            ### ### ### ###
            var_coding_seq = ""
        elif variant_skips_coding_exon:
            ### ### ### ###
            # Simulate coding exon skipping
            # by deletion of coding exon sequence
            ### ### ### ###
            if is_exon_skipped:
                # exon is skipped => delete exon from coding sequence
                logger.debug("Delete the skipped exon from the coding sequence")
                if var_start >= 1:
                    var_coding_seq = (
                        coding_seq[0 : var_start - 1]
                        + coding_seq[var_end : len(coding_seq)]
                    )
                    logger.debug(
                        "Deleting sequence: {}".format(
                            coding_seq[var_start - 1 : var_end]
                        )
                    )
                else:
                    try:
                        assert var_start >= 1
                    except AssertionError:
                        logger.error(
                            "Supplied coding region position should be >= 0\n=> variant position: {}".format(
                                variant_info.to_string()
                            ),
                            exc_info=True,
                        )

                # calculate deletion length and assert its result
                del_length = var_end - var_start + 1
                logger.debug(
                    "len(coding)={}, len(var_coding)={}, del_length={}".format(
                        len(coding_seq), len(var_coding_seq), del_length
                    )
                )
                try:
                    assert len(coding_seq) - len(var_coding_seq) == del_length
                except AssertionError:
                    logger.error(
                        "For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n=> variant position: {}".format(
                            variant_info.to_string()
                        ),
                        exc_info=True,
                    )
            diff_len = -1 * del_length
        else:
            var_coding_seq = coding_seq
    else:
        ### ### ### ###
        # construct observed coding sequence for variant in exonic regions
        ### ### ### ###
        if ">" in var_edit[0:3]:
            logger.debug("Add SNP to coding sequence")
            # get directionally corrected reference and observed bases for SNP
            [ref_seq, obs_seq] = var_edit.split(">")
            if var_start > 1:
                var_coding_seq = (
                    coding_seq[0 : var_start - 1]
                    + obs_seq.lower()
                    + coding_seq[var_start : len(coding_seq)]
                )
            elif var_start == 1:
                var_coding_seq = (
                    obs_seq.lower() + coding_seq[var_start : len(coding_seq)]
                )
            else:
                try:
                    assert var_start >= 0
                except AssertionError:
                    logger.error(
                        "Supplied coding region position should be >= 0\n=> variant position: {}".format(
                            variant_info.to_string()
                        ),
                        exc_info=True,
                    )
            try:
                assert len(coding_seq) == len(var_coding_seq)
            except AssertionError:
                logger.error(
                    "For SNP the sum of length of coding exons should not change\n=> variant position: {}".format(
                        variant_info.to_string()
                    ),
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
            logger.debug("Deleting: {}".format(coding_seq[var_start - 1 : var_end]))
            logger.debug(
                "VEP sequence to delete is: {}".format(transcript_info["var_seq"][0])
            )
            if var_start >= 1:
                var_coding_seq = (
                    coding_seq[0 : var_start - 1]
                    + coding_seq[var_end : len(coding_seq)]
                )
            elif var_start == 0:
                var_coding_seq = coding_seq[var_end : len(coding_seq)]
            else:
                try:
                    assert var_start >= 0
                except AssertionError:
                    logger.error(
                        "Supplied coding region position should be >= 0\n=> variant position: {}".format(
                            variant_info.to_string()
                        ),
                        exc_info=True,
                    )

            # assert deletion by resulted sequence length
            del_length = var_end - var_start + 1
            try:
                logger.debug(
                    "var_coding_seq: {}, len: {}".format(
                        var_coding_seq, len(var_coding_seq)
                    )
                )
                assert len(coding_seq) - len(var_coding_seq) == del_length
            except AssertionError:
                logger.error(
                    "For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n variant position: {}".format(
                        variant_info.to_string()
                    ),
                    exc_info=True,
                )
            ### ### ###
            # introduce the insertion
            ### ### ###
            obs_seq = transcript_info["var_seq"][1]
            try:
                assert len(obs_seq) >= 1
            except AssertionError:
                logger.error(
                    "Supplied VEP does not contain insertion sequence\n=> variant position: {}".format(
                        variant_info.to_string(variant_info.to_string())
                    ),
                    exc_info=True,
                )
            logger.debug("Insert sequence: {}".format(obs_seq))
            if var_start >= 1:
                var_coding_seq_ins = (
                    var_coding_seq[0 : var_start - 1]
                    + obs_seq.lower()
                    + var_coding_seq[var_start - 1 : len(coding_seq)]
                )
            else:
                try:
                    assert var_start >= 1
                except AssertionError:
                    logger.error(
                        "Supplied coding region position should be >= 1\n=> variant position: {}".format(
                            variant_info.to_string()
                        ),
                        exc_info=True,
                    )

            # assert insertion operation by resulted length
            ins_length = len(obs_seq)
            try:
                assert len(var_coding_seq_ins) - len(var_coding_seq) == ins_length
            except AssertionError:
                logger.error(
                    "For insertion should hold: len(var_coding) - len(var_coding with deletion) = len(insertion)\n=> variant position: {}".format(
                        variant_info.to_string()
                    ),
                    exc_info=True,
                )
            var_coding_seq = var_coding_seq_ins
            # assert the length difference for both operations
            diff_len = ins_length - del_length
            try:
                logger.debug(
                    "var_coding_seq: {}, len: {}".format(
                        var_coding_seq, len(var_coding_seq)
                    )
                )
                assert len(var_coding_seq) - len(coding_seq) == diff_len
            except AssertionError:
                logger.error(
                    "For insertion after deletion should hold: len(variant coding) - len(reference coding) = -len(deletion) + len(insertion)\n variant position: {}".format(
                        variant_info.to_string()
                    ),
                    exc_info=True,
                )
        elif "del" in var_edit[0:3]:
            # deletion will be from the start up to end position
            logger.debug("Add deletion to coding sequence")
            if not is_transcript_type_splice_acceptor_donor(
                transcript_info["type_variant"]
            ):
                # assert deleted coding sequence to equal the reference coding sequence
                # include case at which the variant start before or ends after the coding sequence
                if (
                    int(var_coding.pos.start.base) > 0
                    and int(var_coding.pos.end.base) > 0
                ):
                    if not (
                        is_genomic_pos_in_coding_exon(
                            transcript, variant_info.genomic_start
                        )
                        and is_genomic_pos_in_coding_exon(
                            transcript, variant_info.genomic_end
                        )
                    ):
                        # variant start or ends outside an coding exonic range
                        logger.debug(
                            "Variant start or end not in a coding exonic range"
                        )
                        logger.debug(
                            "var start: {}, end: {}".format(var_start, var_end)
                        )
                        affected_coding_length = var_end - var_start + 1
                        logger.debug(
                            "Length of variant sequence that affects coding exonic sequence: {}".format(
                                affected_coding_length
                            )
                        )
                    else:
                        # variant starts and ends inside coding exonic range
                        affected_coding_length = -1
            ### ### ###
            # perform deletion
            ### ### ###
            if var_start >= 1:
                var_coding_seq = (
                    coding_seq[0 : var_start - 1]
                    + coding_seq[var_end : len(coding_seq)]
                )
            else:
                try:
                    assert var_start >= 0
                except AssertionError:
                    logger.error(
                        "Supplied coding region position should be >= 0\n=> variant position: {}".format(
                            variant_info.to_string()
                        ),
                        exc_info=True,
                    )
            del_length = var_end - var_start + 1
            try:
                logger.debug(
                    "var_coding_seq: {}, len: {}".format(
                        var_coding_seq, len(var_coding_seq)
                    )
                )
                assert len(coding_seq) - len(var_coding_seq) == del_length
            except AssertionError:
                logger.error(
                    "For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n variant position: {}".format(
                        variant_info.to_string()
                    ),
                    exc_info=True,
                )
            diff_len = -1 * del_length
        elif "ins" in var_edit[0:3]:
            # for insertion annotation: start & end are the flanking regions
            # the insertion will be placed between the flanking regions
            logger.debug("Add insertion to coding sequence")
            if not (
                is_genomic_pos_in_coding_exon(transcript, variant_info.genomic_start)
                and is_genomic_pos_in_coding_exon(transcript, variant_info.genomic_end)
            ):
                # variant start or ends outside an coding exonic range
                logger.debug("Variant start or end not in a coding exonic range")
                affected_coding_length = var_end - var_start + 1
                logger.debug(
                    "Length of variant sequence that affects coding exonic sequence: {}".format(
                        affected_coding_length
                    )
                )
            else:
                # variant inside coding exonic range
                affected_coding_length = -1

            ### ### ###
            # calculate sequence to be inserted
            ### ### ###
            if len(transcript_info["var_seq"][0]) > 0:
                # VEP contains insertion sequence
                if "_" in transcript_info["var_seq"][0]:
                    # insertion sequence is described by coding coordinates
                    logger.info(
                        "Insertion sequence is described by coding coordinates \n=> variant pos: {}".format(
                            variant_info.to_string()
                        )
                    )
                    [ins_source_start, ins_source_end] = transcript_info["var_seq"][
                        0
                    ].split("_")
                    transcript_info["var_seq"][0] = coding_seq[
                        int(ins_source_start) - 1 : int(ins_source_end.strip())
                    ]
                if affected_coding_length == -1:
                    # all insertion inside coding exon
                    obs_seq = transcript_info["var_seq"][0]
                else:
                    obs_seq = transcript_info["var_seq"][0][-affected_coding_length:]
            else:
                try:
                    assert 1 == 0
                except AssertionError:
                    logger.error(
                        "vep does not contain insertion sequence\n=> variant position: {}".format(
                            variant_info.to_string()
                        ),
                        exc_info=True,
                    )
            logger.debug("Insert the sequence: {}".format(obs_seq))
            ### ### ###
            # perform insertion
            ### ### ###
            if var_start >= 1:
                if var_start == var_end:
                    # start equal ends so the rest part, after the insertion, should start on end position
                    if vep_coordinates_used:
                        var_coding_seq = (
                            coding_seq[0:var_start]
                            + obs_seq.lower()
                            + coding_seq[var_end : len(coding_seq)]
                        )
                    else:
                        var_coding_seq = (
                            coding_seq[0:var_start]
                            + obs_seq.lower()
                            + coding_seq[var_end : len(coding_seq)]
                        )
                else:
                    if vep_coordinates_used:
                        var_coding_seq = (
                            coding_seq[0:var_start]
                            + obs_seq.lower()
                            + coding_seq[var_end - 1 : len(coding_seq)]
                        )
                    else:
                        var_coding_seq = (
                            coding_seq[0:var_start]
                            + obs_seq.lower()
                            + coding_seq[var_end - 1 : len(coding_seq)]
                        )
            else:
                try:
                    assert var_start >= 1
                except AssertionError:
                    logger.error(
                        "Supplied coding region position should be >= 0\n=> variant position: {}".format(
                            variant_info.to_string()
                        ),
                        exc_info=True,
                    )
            ins_length = len(obs_seq)
            try:
                assert len(var_coding_seq) - len(coding_seq) == ins_length
            except AssertionError:
                logger.error(
                    "For insertion should hold: len(var_coding) - len(reference_coding) = len(insertion)\n=> variant position: {}".format(
                        variant_info.to_string()
                    ),
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
                is_genomic_pos_in_coding_exon(transcript, variant_info.genomic_start)
                and is_genomic_pos_in_coding_exon(transcript, variant_info.genomic_end)
            ):
                ### ### ###
                # variant starts or ends in intronic region
                ### ### ###
                logger.debug("Variant start or end not in a coding exonic range")
                logger.debug("var start: {}, end: {}".format(var_start, var_end))
                affected_coding_length = var_end - var_start + 1
                logger.debug(
                    "Length of variant sequence that affects coding exonic sequence: {}".format(
                        affected_coding_length
                    )
                )
            else:
                ### ### ###
                # variant inside coding exonic range
                ### ### ###
                affected_coding_length = -1

            ### ### ###
            # calculate sequence to be duplicated
            ### ### ###
            if len(transcript_info["var_seq"][0]) > 0:
                if affected_coding_length == -1:
                    # all duplication inside coding exon
                    obs_seq = transcript_info["var_seq"][0]
                else:
                    obs_seq = transcript_info["var_seq"][0][-affected_coding_length:]
            else:
                if var_start == var_end:
                    obs_seq = coding_seq[var_start - 1]
                else:
                    obs_seq = coding_seq[var_start - 1 : var_end]
            logger.debug("Duplicate the sequence: {}".format(obs_seq))

            ### ### ###
            # perform duplication
            ### ### ###
            if var_start == var_end:
                # duplication of one nucleotide
                logger.debug("Duplication of one nucleotide")
                if var_start >= 1:
                    var_coding_seq = (
                        coding_seq[0:var_start]
                        + obs_seq.lower()
                        + coding_seq[var_start : len(coding_seq)]
                    )
                else:
                    try:
                        assert var_start >= 0
                    except AssertionError:
                        logger.error(
                            "Supplied coding region position should be >= 0\n=> variant position: {}".format(
                                variant_info.to_string()
                            ),
                            exc_info=True,
                        )
            else:
                # duplication of multiple nucleotides
                logger.debug("Duplication of multiple nucleotides")
                if var_start >= 1:
                    if vep_coordinates_used:
                        var_coding_seq = (
                            coding_seq[0:var_end]
                            + obs_seq.lower()
                            + coding_seq[var_end : len(coding_seq)]
                        )
                    else:
                        var_coding_seq = (
                            coding_seq[0 : var_end + 1]
                            + obs_seq.lower()
                            + coding_seq[var_end + 1 : len(coding_seq)]
                        )
                else:
                    try:
                        assert var_start >= 1
                    except AssertionError:
                        logger.error(
                            "Supplied coding region position should be >= 0\n=> variant position: {}".format(
                                variant_info.to_string()
                            ),
                            exc_info=True,
                        )
            dupl_length = len(obs_seq)
            try:
                assert len(var_coding_seq) - len(coding_seq) == dupl_length
            except AssertionError:
                logger.error(
                    "For duplication should hold: len(var_coding) - len(reference_coding) = len(duplication)\n variant position: {}".format(
                        variant_info.to_string()
                    ),
                    exc_info=True,
                )
                diff_len = +1 * dupl_length

    ### ### ### ###
    # after creating coding sequence with variant,
    # check that indeed is different from the reference
    ### ### ### ###
    if variant_skips_coding_exon or not is_transcript_type_splice_acceptor_donor(
        transcript_info["type_variant"]
    ):
        try:
            assert coding_seq.upper() != var_coding_seq.upper()
        except AssertionError:
            logger.error(
                "Coding sequence of reference and sample should be different\n variant position: {}".format(
                    variant_info.to_string()
                ),
                exc_info=True,
            )
        print_ref_observed_seq(
            coding_seq,
            var_coding_seq,
            transcript_info,
            var_start,
            var_end,
            var_edit,
            is_exon_skipped,
        )

        ### ### ### ### ### ###
        # assert constructed variant coding sequence
        # contains start (if edit not start_lost) and stop codon
        ### ### ### ### ### ###
        if variant_skips_start_codon_exon and variant_skips_stop_codon_exon:
            pass
        elif variant_skips_start_codon_exon:
            try:
                assert var_coding_seq[-3:] in ["TAG", "TAA", "TGA"]
            except AssertionError:
                logger.error(
                    "Constructed observed coding sequence should contain stop codon\n=> variant position: {}".format(
                        variant_info.to_string()
                    ),
                    exc_info=True,
                )
        elif variant_skips_stop_codon_exon:
            if coding_seq[0:3] == "ATG":
                try:
                    assert var_coding_seq[0:3] == "ATG"
                except AssertionError:
                    logger.debug(
                        "Constructed observed coding sequence should contain start codon\n=> variant position: {}".format(
                            variant_info.to_string()
                        )
                    )
        else:
            ### ### ###
            # variant does not skip exon => prepare range of variant start stop
            # to assert for existence of start and stop codons
            ### ### ###
            if var_start == var_end:
                var_positions = range(var_start, var_start + 1)
            else:
                var_positions = range(var_start, var_end)

            if not (
                "start_lost" in transcript_info["type_variant"]
                or "start_retained" in variant_info.variant_type
            ):
                start_positions = set(range(1, 4))
                if (
                    not start_positions.intersection(var_positions)
                    and coding_seq[0:3] == "ATG"
                ):
                    # if variant does not change positions on the start of the coding sequence
                    # and ensembl record for transcript starts by ATG
                    # assert that the observed coding sequence starts with ATG
                    try:
                        assert var_coding_seq[0:3] == "ATG"
                    except AssertionError:
                        logger.error(
                            "Constructed observed coding sequence should contain start codon\n=> variant position: {}".format(
                                variant_info.to_string()
                            ),
                            exc_info=True,
                        )

                if not (
                    "stop_retained_variant" in transcript_info["type_variant"]
                    or "stop_lost" in transcript_info["type_variant"]
                ):
                    stop_positions = set(
                        range(len(coding_seq) + 1 - 3, len(coding_seq) + 1)
                    )
                    if not stop_positions.intersection(var_positions):
                        # if variant does not change positions in the end of the coding sequence
                        # assert that the observed coding sequence finishes with "TAG", "TAA" or "TGA"
                        try:
                            assert var_coding_seq[-3:] in ["TAG", "TAA", "TGA"]
                        except AssertionError:
                            logger.error(
                                "Constructed observed coding sequence should contain stop codon\n=> variant position: {}".format(
                                    variant_info.to_string()
                                ),
                                exc_info=True,
                            )

    return var_coding_seq.upper(), diff_len


def print_ref_observed_seq(
    coding_seq,
    var_coding_seq,
    transcript_info,
    var_start,
    var_end,
    var_edit,
    is_exon_skipped,
):
    """
    Print reference and observed coding sequence

    Parameters
    ----------
    coding_seq : str
        reference coding sequence
    var_coding_seq : str
        observed coding sequence
    transcript_info : dict of str: str
        variant-affected transcript GSvar information
    var_start : int
        variant start position (in coding sequence)
    var_end : int
        variant end position (in coding sequence)
    var_edit : str
        variant edit information
    is_exon_skipped : bool
        variant causes an exon to be skipped (True), otherwise (False)

    Returns
    -------
    None
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
        is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"])
        and is_exon_skipped
    ):
        # deletion or skipped exon case
        ref_deleted_var = (
            coding_seq[var_start - print_start - 1 : var_start - 1]
            + coding_seq[var_start - 1 : var_end].lower()
            + coding_seq[var_end : var_end + print_end]
        )
        logger.debug(">coding seq:\n {}".format(ref_deleted_var))
    else:
        # insertion, duplication, SNP
        logger.debug(
            ">coding seq:\n {}".format(
                coding_seq[var_start - print_start - 1 : var_end + print_end]
            )
        )

    if var_start == 0:
        logger.debug(
            ">var coding seq:\n {}".format(
                var_coding_seq[var_start : var_end + print_end]
            )
        )
    else:
        logger.debug(
            ">var coding seq:\n {}".format(
                var_coding_seq[var_start - print_start - 1 : var_end + print_end]
            )
        )
