#!/usr/bin/env python3

from refactoring.variant import VariantInfo, TranscriptInfo

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_example_splice_acceptor_BRCA1() -> tuple:
    """
    Create BRCA1 splice acceptor variant
    On - strand, in splice_site
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5278-1G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        ["splice_acceptor"],
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


def create_example_splice_donor_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    On - strand, in splice_site
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5332+2T>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        ["splice_donor"],
        hgvs,
        5332,
        5332,
        exon=None,
        intron=20,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        ["splice_donor"],
        "17",
        41203078,
        41203078,
        "some_id",
        "A",
        "G",
    )
    return (transcript, variant)


def create_example_splice_donor_exonic_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    On - strand, in splice_site
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4484G>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        ["splice_donor"],
        hgvs,
        4484,
        4484,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        ["splice_donor"],
        "17",
        41228505,
        41228505,
        "some_id",
        "C",
        "G",
    )
    return (transcript, variant)


def create_example_splice_acceptor_exonic_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    On - strand, in splice_site
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4484G>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        ["splice_donor"],
        hgvs,
        4484,
        4484,
        exon=None,
        intron=20,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        ["splice_donor"],
        "17",
        41228505,
        41228505,
        "some_id",
        "C",
        "G",
    )
    return (transcript, variant)


def create_example_splice_donor_BRCA1_2() -> tuple:
    """
    Create BRCA1 splice donor variant
    From ClinVar: VCV125738.12
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4986+1G>T".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        ["splice_donor"],
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


def create_example_splice_BARD1() -> tuple:
    """
    Create BARD1 splice donor variant
    From ClinVar: VCV125738.12
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1569-7T>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000260947",
        [
            "splice_region_variant",
            "splice_polypyrimidine_tract_variant",
            "intron_variant",
        ],
        hgvs,
        1569,
        1569,
        exon=None,
        intron=6,
        var_protein=None,
    )
    variant = VariantInfo(
        "BARD1",
        ["splice_donor"],
        "2",
        215617286,
        215617286,
        "some_id",
        "A",
        "C",
    )
    return (transcript, variant)


def create_example_splice_donor_TP53() -> tuple:
    """
    Create TP53 splice donor variant
    On - strand, outside splice_site
    From ClinVar: VCV125738.12
    """
    hgvs = hgvs_parser.parse_c_posedit("c.357+5G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000269305",
        ["intron_variant", "splice_donor_5th_base_variant"],
        hgvs,
        375,
        375,
        exon=None,
        intron=4,
        var_protein=None,
    )
    variant = VariantInfo(
        "TP53",
        ["intron_variant", "splice_donor_5th_base_variant"],
        "17",
        7579307,
        7579307,
        "some_id",
        "C",
        "T",
    )
    return (transcript, variant)


def create_example_splice_acceptor_BRCA2() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000038132.31
    """
    hgvs = hgvs_parser.parse_c_posedit("c.7977-1G>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        ["splice_acceptor"],
        hgvs,
        7977,
        7977,
        exon=None,
        intron=17,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        ["splice_acceptor"],
        "13",
        32937315,
        32937315,
        "some_id",
        "G",
        "C",
    )
    return (transcript, variant)


def create_example_splice_acceptor_exonic_BRCA2() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000052458.15
    """
    hgvs = hgvs_parser.parse_c_posedit("c.7978T>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        ["splice_acceptor", "missense_variant"],
        hgvs,
        7978,
        7978,
        exon=18,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        ["splice_acceptor", "missense_variant"],
        "13",
        32937317,
        32937317,
        "some_id",
        "T",
        "G",
    )
    return (transcript, variant)


def create_example_splice_donor_BRCA2() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000267692.60
    """
    hgvs = hgvs_parser.parse_c_posedit("c.8331+2T>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        ["splice_donor"],
        hgvs,
        8331,
        8331,
        exon=None,
        intron=18,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        ["splice_donor"],
        "13",
        32937672,
        32937672,
        "some_id",
        "T",
        "C",
    )
    return (transcript, variant)


def create_example_splice_donor_exonic_BRCA2() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000052552.2
    """
    hgvs = hgvs_parser.parse_c_posedit("c.8331G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        ["splice_donor", "synonymous_variant"],
        hgvs,
        8331,
        8331,
        exon=18,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        ["splice_donor", "synonymous_variant"],
        "13",
        32937670,
        32937670,
        "some_id",
        "G",
        "A",
    )
    return (transcript, variant)


def create_example_splice_variant_BRCA2() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000220547.6
    """
    hgvs = hgvs_parser.parse_c_posedit("c.7618-10T>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        ["splice_region_variant", "intron_variant"],
        hgvs,
        7618,
        7618,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        ["splice_region_variant", "intron_variant"],
        "13",
        32931869,
        32931869,
        "some_id",
        "T",
        "G",
    )
    return (transcript, variant)


def create_example_splice_variant_BRCA2_2() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000918812.2
    """
    hgvs = hgvs_parser.parse_c_posedit("c.7805+8A>c".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        ["splice_region_variant", "intron_variant"],
        hgvs,
        7805,
        7805,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        ["splice_region_variant", "intron_variant"],
        "13",
        32932074,
        32932074,
        "some_id",
        "A",
        "C",
    )
    return (transcript, variant)
