#!/usr/bin/env python3

from refactoring.variant import (
    PopulationDatabases_gnomAD,
    VariantInfo,
    TranscriptInfo,
    Variant,
    PopulationDatabases,
    AffectedRegion,
    VARTYPE,
)

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_test_variant() -> Variant:
    transinfo, varinfo = create_example_splice_acceptor_BRCA1()
    region = AffectedRegion(True, False)
    flossies = PopulationDatabases("Flossies", 0.5)
    gnomad = PopulationDatabases_gnomAD("gnomAd", 0.05, 20, "AES", 0.06, 10)
    prediction = {"SpliceAI": 0.5, "REVEL": 0.5}
    return Variant(varinfo, [transinfo], prediction, gnomad, flossies, region)


def create_test_variant_no_flossies() -> Variant:
    transinfo, varinfo = create_example_splice_acceptor_BRCA1()
    region = AffectedRegion(True, False)
    gnomad = PopulationDatabases_gnomAD("gnomAd", 0.05, 20, "AES", 0.06, 10)
    prediction = {"SpliceAI": 0.5, "REVEL": 0.5}
    return Variant(
        variant_info=varinfo,
        transcript_info=[transinfo],
        prediction_tools=prediction,
        gnomad=gnomad,
        affected_region=region,
    )


def create_test_variant_no_hotspot() -> Variant:
    transinfo, varinfo = create_example_splice_acceptor_BRCA1()
    region = AffectedRegion(repetitive_region=True)
    gnomad = PopulationDatabases_gnomAD("gnomAd", 0.05, 20, "AES", 0.06, 10)
    flossies = PopulationDatabases("Flossies", 0.5)
    prediction = {"SpliceAI": 0.5, "REVEL": 0.5}
    return Variant(
        variant_info=varinfo,
        transcript_info=[transinfo],
        prediction_tools=prediction,
        gnomad=gnomad,
        flossies=flossies,
        affected_region=region,
    )


def create_example_splice_acceptor_BRCA1() -> tuple:
    """
    Create BRCA1 splice acceptor variant
    On - strand, in splice_site
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5278-1G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.SPLICE_ACCEPTOR],
        hgvs,
        5278,
        5278,
        exon=None,
        intron=19,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        [VARTYPE.SPLICE_ACCEPTOR],
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
        [VARTYPE.SPLICE_DONOR],
        hgvs,
        5332,
        5332,
        exon=None,
        intron=20,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
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
        [VARTYPE.SPLICE_DONOR],
        hgvs,
        4484,
        4484,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
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
        [VARTYPE.SPLICE_DONOR],
        hgvs,
        4484,
        4484,
        exon=None,
        intron=20,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
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
        [VARTYPE.SPLICE_DONOR],
        hgvs,
        4986,
        4986,
        exon=None,
        intron=15,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
        "17",
        41222944,
        41222944,
        "some_id",
        "C",
        "T",
    )
    return (transcript, variant)


def create_example_minus_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    - strand, outside of splice site
    From ClinVar: VCV000433731.4
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5278-6T>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.SPLICE_DONOR],
        hgvs,
        5278,
        5278,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
        "17",
        41203140,
        41203140,
        "some_id",
        "A",
        "G",
    )
    return (transcript, variant)


def create_example_plus_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    - strand, outside of splice site
    From ClinVar: VCV000125827.18
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5332+4A>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.SPLICE_DONOR],
        hgvs,
        5332,
        5332,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
        "17",
        41203140,
        41203140,
        "some_id",
        "T",
        "C",
    )
    return (transcript, variant)


def create_example_splice_BARD1() -> tuple:
    """
    Create BARD1 splice donor variant
    - strand,
    From ClinVar: VCV125738.12
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1569-7T>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000260947",
        [
            VARTYPE.SPLICE_REGION_VARIANT,
            VARTYPE.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
            VARTYPE.INTRON_VARIANT,
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
        [
            VARTYPE.SPLICE_REGION_VARIANT,
            VARTYPE.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
            VARTYPE.INTRON_VARIANT,
        ],
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
        [VARTYPE.INTRON_VARIANT, VARTYPE.SPLICE_DONOR_5TH_BASE_VARIANT],
        hgvs,
        375,
        375,
        exon=None,
        intron=4,
        var_protein=None,
    )
    variant = VariantInfo(
        "TP53",
        [VARTYPE.INTRON_VARIANT, VARTYPE.SPLICE_DONOR_5TH_BASE_VARIANT],
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
        [VARTYPE.SPLICE_ACCEPTOR],
        hgvs,
        7977,
        7977,
        exon=None,
        intron=17,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        [VARTYPE.SPLICE_ACCEPTOR],
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
        [VARTYPE.SPLICE_ACCEPTOR, VARTYPE.MISSENSE_VARIANT],
        hgvs,
        7978,
        7978,
        exon=18,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        [VARTYPE.SPLICE_ACCEPTOR, VARTYPE.MISSENSE_VARIANT],
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
        [VARTYPE.SPLICE_DONOR],
        hgvs,
        8331,
        8331,
        exon=None,
        intron=18,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        [VARTYPE.SPLICE_DONOR],
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
        [VARTYPE.SPLICE_DONOR, VARTYPE.SYNONYMOUS_VARIANT],
        hgvs,
        8331,
        8331,
        exon=18,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        [VARTYPE.SPLICE_DONOR, VARTYPE.SYNONYMOUS_VARIANT],
        "13",
        32937670,
        32937670,
        "some_id",
        "G",
        "A",
    )
    return (transcript, variant)


def create_example_splice_variant_BRCA2_minus() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000220547.6
    """
    hgvs = hgvs_parser.parse_c_posedit("c.7618-10T>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        [VARTYPE.SPLICE_REGION_VARIANT, VARTYPE.INTRON_VARIANT],
        hgvs,
        7618,
        7618,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        [VARTYPE.SPLICE_REGION_VARIANT, VARTYPE.INTRON_VARIANT],
        "13",
        32931869,
        32931869,
        "some_id",
        "T",
        "G",
    )
    return (transcript, variant)


def create_example_splice_variant_BRCA2_plus() -> tuple:
    """
    Create BRCA2 splice acceptor variant
    On + strand, in splice_site
    From ClinVar: VCV000918812.2
    """
    hgvs = hgvs_parser.parse_c_posedit("c.7805+8A>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        [VARTYPE.SPLICE_REGION_VARIANT, VARTYPE.INTRON_VARIANT],
        hgvs,
        7805,
        7805,
        exon=None,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "BRCA2",
        [VARTYPE.SPLICE_REGION_VARIANT, VARTYPE.INTRON_VARIANT],
        "13",
        32932074,
        32932074,
        "some_id",
        "A",
        "C",
    )
    return (transcript, variant)
