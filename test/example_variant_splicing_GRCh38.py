#!/usr/bin/env python3

from variant_classification.variant import (
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
    flossies = PopulationDatabases("Flossies", count=5, frequency=None)
    gnomad = PopulationDatabases_gnomAD("gnomAd", 0.05, 20, "AES", 0.06, 10)
    prediction = {"SpliceAI": 0.5, "REVEL": 0.5}
    return Variant(
        variant_info=varinfo,
        transcript_info=[transinfo],
        prediction_tools=prediction,
        gnomad=gnomad,
        flossies=flossies,
        affected_region=region,
    )


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
    region = AffectedRegion(critical_region=True, cancer_hotspot=None, cold_spot=None)
    gnomad = PopulationDatabases_gnomAD("gnomAd", 0.05, 20, "AES", 0.06, 10)
    flossies = PopulationDatabases("Flossies", 0.5, None)
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
    Classified as pathogenic (PVS1, PS1, PM2)
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
        "17",
        43051118,
        43051118,
        "BRCA1",
        [VARTYPE.SPLICE_ACCEPTOR],
        "C",
        "T",
    )
    return (transcript, variant)


def create_example_splice_donor_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    On - strand, in splice_site
    From HerediVar
    No classification information
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5332+1G>C".split("c.")[1])
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
        "17",
        43051062,
        43051062,
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
        "C",
        "G",
    )
    return (transcript, variant)


def create_example_splice_donor_exonic_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    On - strand, in splice_site
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4484G>T".split("c.")[1])
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
        "17",
        43076488,
        43076488,
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
        "C",
        "A",
    )
    return (transcript, variant)


def create_example_splice_acceptor_exonic_BRCA1() -> tuple:
    """
    Create BRCA1 splice donor variant
    On - strand, in splice_site
    From HerediVar
    Not classified
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
        "17",
        41228505,
        41228505,
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
        "C",
        "G",
    )
    return (transcript, variant)


def create_example_splice_donor_BRCA1_2() -> tuple:
    """
    Create BRCA1 splice donor variant
    Previously classified in HerediVar as 4
    From ClinVar: VCV000125738.12
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4986+1G>A".split("c.")[1])
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
        "17",
        43070927,
        43070927,
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
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
        "17",
        43051123,
        43051123,
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
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
        "17",
        43051059,
        43051059,
        "BRCA1",
        [VARTYPE.SPLICE_DONOR],
        "T",
        "C",
    )
    return (transcript, variant)


def create_example_splice_BARD1() -> tuple:
    """
    Create BARD1 splice donor variant
    - strand,
    From ClinVar: VCV001662301.5
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1569-7T>C".split("c.")[1])
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
        "2",
        214752562,
        214752562,
        "BARD1",
        [
            VARTYPE.SPLICE_REGION_VARIANT,
            VARTYPE.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
            VARTYPE.INTRON_VARIANT,
        ],
        "A",
        "G",
    )
    return (transcript, variant)


def create_example_splice_donor_TP53() -> tuple:
    """
    Create TP53 splice donor variant
    On - strand, outside splice_site
    From ClinVar: VCV000481015.14
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
        "17",
        7675989,
        7675989,
        "TP53",
        [VARTYPE.INTRON_VARIANT, VARTYPE.SPLICE_DONOR_5TH_BASE_VARIANT],
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
        "13",
        32363178,
        32363178,
        "BRCA2",
        [VARTYPE.SPLICE_ACCEPTOR],
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
        "13",
        32363180,
        32363180,
        "BRCA2",
        [VARTYPE.SPLICE_ACCEPTOR, VARTYPE.MISSENSE_VARIANT],
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
        "13",
        32363535,
        32363535,
        "BRCA2",
        [VARTYPE.SPLICE_DONOR],
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
        "13",
        32363533,
        32363533,
        "BRCA2",
        [VARTYPE.SPLICE_DONOR, VARTYPE.SYNONYMOUS_VARIANT],
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
        "13",
        32357735,
        32357735,
        "BRCA2",
        [VARTYPE.SPLICE_REGION_VARIANT, VARTYPE.INTRON_VARIANT],
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
        "13",
        32357937,
        32357937,
        "BRCA2",
        [VARTYPE.SPLICE_REGION_VARIANT, VARTYPE.INTRON_VARIANT],
        "A",
        "C",
    )
    return (transcript, variant)
