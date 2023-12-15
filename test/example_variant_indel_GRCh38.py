#!/usr/bin/env python3

from variant_classification.variant import VariantInfo, TranscriptInfo, VARTYPE

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_example_dup() -> tuple:
    """
    Create PMS2 duplication variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1831dup".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000265849",
        [VARTYPE.FRAMESHIFT_VARIANT],
        hgvs,
        1831,
        1832,
        exon=11,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "7",
        5986933,
        5986934,
        "PMS2",
        [VARTYPE.FRAMESHIFT_VARIANT],
        "A",
        "AT",
    )
    return (transcript, variant)


def create_example_ins() -> tuple:
    """
    Create BRCA1 insertion variant
    From ClinVar: VCV000266352.6
    """
    hgvs = hgvs_parser.parse_c_posedit("c.330_331insA".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.FRAMESHIFT_VARIANT],
        hgvs,
        330,
        331,
        exon=6,
        intron=None,
        var_protein="p.Glu111fs",
    )
    variant = VariantInfo(
        "17",
        43104232,
        43104233,
        "BRCA1",
        [VARTYPE.FRAMESHIFT_VARIANT],
        "",
        "T",
    )
    return (transcript, variant)


def create_example_ins_no_frameshift() -> tuple:
    """
    Create BRCA1 insertion variant that leads to no frameshift
    Altered from ClinVar: VCV000266352.6
    """
    hgvs = hgvs_parser.parse_c_posedit("c.330_331insAAA".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.FRAMESHIFT_VARIANT],
        hgvs,
        330,
        331,
        exon=6,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "17",
        43104232,
        43104233,
        "BRCA1",
        [VARTYPE.FRAMESHIFT_VARIANT],
        "",
        "TTT",
    )
    return (transcript, variant)


def create_example_del() -> tuple:
    """
    Create TP53 deletion variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.286_288del".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000269305",
        [VARTYPE.INFRAME_DELETION],
        hgvs,
        286,
        288,
        exon=4,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "17",
        7676080,
        7676082,
        "TP53",
        [VARTYPE.INFRAME_DELETION],
        "CAGA",
        "C",
    )
    return (transcript, variant)


def create_example_del_frameshift() -> tuple:
    """
    Create CDH1 deletion variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1567del".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000404276",
        [VARTYPE.INFRAME_DELETION],
        hgvs,
        1567,
        1568,
        exon=15,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "22",
        28687961,
        28687962,
        "CHEK2",
        [VARTYPE.FRAMESHIFT_VARIANT],
        "G",
        "",
    )
    return (transcript, variant)


def create_example_indel() -> tuple:
    """
    Create BRCA1 indel variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4391_4393delinsTT".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.FRAMESHIFT_VARIANT],
        hgvs,
        4391,
        4393,
        exon=13,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "17",
        43076579,
        43076581,
        "BRCA1",
        [VARTYPE.FRAMESHIFT_VARIANT],
        "some_id",
        "TAG",
    )
    return (transcript, variant)
