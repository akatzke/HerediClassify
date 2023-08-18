#!/usr/bin/env python3

from refactoring.variant import VariantInfo, TranscriptInfo

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_example_splice_acceptor() -> tuple:
    """
    Create BRCA1 splice acceptor variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5278-1G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        "splice_acceptor",
        hgvs,
        5278,
        5278,
        exon=None,
        intron=19,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        ["splice_acceptor"],
        "17",
        41203135,
        41203135,
        "some_id",
        "C",
        "T",
    )
    return (transcript, variant)


def create_example_splice_donor() -> tuple:
    """
    Create BRCA1 splice donor variant
    From ClinVar: VCV125738.12
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4986+1G>T".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        "splice_donor",
        hgvs,
        4986,
        4986,
        exon=None,
        intron=15,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        ["splice_donor"],
        "17",
        41222944,
        41222944,
        "some_id",
        "C",
        "T",
    )
    return (transcript, variant)


def create_example_dup() -> tuple:
    """
    Create PMS2 duplication variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1831dup".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000265849",
        "frameshift_variant",
        hgvs,
        1831,
        1832,
        exon=11,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "PMS2",
        ["frameshift_variant"],
        "7",
        6026564,
        6026565,
        "some_id",
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
        "frameshift_variant",
        hgvs,
        330,
        331,
        exon=6,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1", ["frameshift_variant"], "17", 41256249, 41256250, "some_id", "", "T"
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
        "frameshift_variant",
        hgvs,
        330,
        331,
        exon=6,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1", ["frameshift_variant"], "17", 41256249, 41256250, "some_id", "", "TTT"
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
        "inframe_deletion",
        hgvs,
        286,
        288,
        exon=4,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "TP53", ["inframe_deletion"], "17", 7579399, 7579401, "some_id", "CAGA", "C"
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
        "frameshift_variant",
        hgvs,
        4391,
        4393,
        exon=13,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        ["frameshift_variant"],
        "17",
        41228596,
        41228598,
        "some_id",
        "TAG",
        "AA",
    )
    return (transcript, variant)


def create_example_snv_mis() -> tuple:
    """
    Create BRCA1 missense SNV variant
    From HerediVar
    """
    hgvs_1 = hgvs_parser.parse_c_posedit("c.5219T>G".split("c.")[1])
    transcript_1 = TranscriptInfo(
        "ENST00000357654",
        "missense_variant",
        hgvs_1,
        5219,
        5219,
        exon=19,
        intron=None,
        var_protein="p.Val1740Gly",
    )
    hgvs_2 = hgvs_parser.parse_c_posedit("c.5282T>G".split("c.")[1])
    transcript_2 = TranscriptInfo(
        "ENST00000471181",
        "missense_variant",
        hgvs_2,
        5282,
        5282,
        exon=20,
        intron=None,
        var_protein="p.Val1761Gly",
    )
    transcripts = [transcript_1, transcript_2]
    variant = VariantInfo(
        "BRCA1", ["missense_variant"], "17", 41209127, 41209127, "some_id", "A", "C"
    )
    return (transcripts, variant)


def create_example_snv_syn() -> tuple:
    """
    Create BRCA1 synonymous SNV variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1740C>T".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        "synonymous_variant",
        hgvs,
        4986,
        4986,
        exon=10,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1", ["synonymous_variant"], "17", 41245808, 4125808, "some_id", "G", "A"
    )
    return (transcript, variant)


def create_start_lost() -> tuple:
    """
    Create RAD51D start lost variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.2T>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000345365",
        "start lost",
        hgvs,
        2,
        2,
        exon=1,
        intron=None,
        var_protein="p.Met1?",
    )
    variant = VariantInfo(
        "RAD51D", ["start_lost"], "17", 33446631, 33446631, "some_id", "A", "C"
    )
    return (transcript, variant)


def create_brca2_missense() -> tuple:
    """
    Variant defined in origianl GenOtoScope test set
    """
    hgvs = hgvs_parser.parse_c_posedit("c.3073A>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        "missense_variant",
        hgvs,
        3073,
        3073,
        exon=11,
        intron=None,
        var_protein="p.Lys1025Glu",
    )
    variant = VariantInfo(
        "BRCA2",
        ["missense_variant"],
        "13",
        32911565,
        32911565,
        "RS=rs80358550",
        "A",
        "G",
    )
    return (transcript, variant)


def create_myo6_frameshift() -> tuple:
    """
    Variant defined in origianl GenOtoScope test set
    """
    hgvs = hgvs_parser.parse_c_posedit("c.2854dup".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000369975",
        "frameshift_variant",
        hgvs,
        2854,
        2854,
        exon=25,
        intron=None,
        var_protein="p.Glu952GlyfsTer8",
    )
    variant = VariantInfo(
        "MYO6",
        ["frameshift_variant"],
        "6",
        76599967,
        76599967,
        "some_id",
        "A",
        "AG",
    )
    return (transcript, variant)


def create_diaph1_frameshift() -> tuple:
    """
    Variant defined in origianl GenOtoScope test set
    """
    hgvs = hgvs_parser.parse_c_posedit("c.3627_3628del".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000253811",
        "frameshift_variant",
        hgvs,
        3627,
        3628,
        exon=27,
        intron=None,
        var_protein="p.Ala1211SerfsTer31",
    )
    variant = VariantInfo(
        "DIAPH1",
        ["frameshift_variant"],
        "5",
        140903745,
        140903745,
        "some_id",
        "CCT",
        "C",
    )
    return (transcript, variant)
