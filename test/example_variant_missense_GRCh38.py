#!/usr/bin/env python3

from refactoring.variant import VariantInfo, TranscriptInfo, VARTYPE

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_example_mis() -> tuple:
    """
    Create BRCA1 missense SNV variant
    - strand, 2 position in codon, VUS no matching ClinVar
    From HerediVar
    """
    hgvs_1 = hgvs_parser.parse_c_posedit("c.5219T>G".split("c.")[1])
    transcript_1 = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.MISSENSE_VARIANT],
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
        [VARTYPE.MISSENSE_VARIANT],
        hgvs_2,
        5282,
        5282,
        exon=20,
        intron=None,
        var_protein="p.Val1761Gly",
    )
    transcripts = [transcript_1, transcript_2]
    variant = VariantInfo(
        "17",
        43057110,
        43057110,
        "BRCA1",
        [VARTYPE.MISSENSE_VARIANT],
        "A",
        "C",
    )
    return (transcripts, variant)


def create_example_mis_first_base() -> tuple:
    """
    Create BRCA1 missense SNV variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5218G>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        5218,
        5218,
        exon=19,
        intron=None,
        var_protein="p.Val1740Leu",
    )
    variant = VariantInfo(
        "17",
        43057111,
        43057111,
        "BRCA1",
        [VARTYPE.MISSENSE_VARIANT],
        "C",
        "G",
    )
    return (transcript, variant)


def create_example_mis_BRCA1_2() -> tuple:
    """
    Create BRCA1 missense SNV variant
    - strand, 2 position in codon, matching pathogenic ClinVar entries
    From HerediVar
    Classified as class (PS3, PM1, PM2_pp, PP1, BP4)
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4232T>C".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        4232,
        4232,
        exon=12,
        intron=None,
        var_protein="p.Met1411Thr",
    )
    variant = VariantInfo(
        "17",
        43082529,
        43082529,
        "BRCA1",
        [VARTYPE.MISSENSE_VARIANT],
        "A",
        "G",
    )
    return (transcript, variant)


def create_example_mis_BRCA1_1() -> tuple:
    """
    Create BRCA1 missense SNV variant
    - strand, 1 position in codon, matching pathogenic ClinVar entries
    From ClinVar: VCV000141465.12
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4231A>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        4231,
        4231,
        exon=12,
        intron=None,
        var_protein="p.Met1411Thr",
    )
    variant = VariantInfo(
        "17",
        43082530,
        43082530,
        "BRCA1",
        [VARTYPE.MISSENSE_VARIANT],
        "T",
        "C",
    )
    return (transcript, variant)


def create_example_mis_BRCA1_3() -> tuple:
    """
    Create BRCA1 missense SNV variant
    - strand, 3 position in codon, matching pathogenic ClinVar entries
    From ClinVar: VCV000628271.4
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4233G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        4233,
        4233,
        exon=12,
        intron=None,
        var_protein="p.Met1411Ile",
    )
    variant = VariantInfo(
        "17",
        43082528,
        43082528,
        "BRCA1",
        [VARTYPE.MISSENSE_VARIANT],
        "C",
        "T",
    )
    return (transcript, variant)


def create_example_mis_BRCA1_2_exons() -> tuple:
    """
    Create BRCA1 missense SNV variant
    - strand, 3 position in codon, matching pathogenic ClinVar entries
    Codon spans multiple exons
    Codon genomic pos: [41243048, 41243049, 41243452]
    ClinVar ids: [VCV001737791.3, VCV000234236.13, VCV000096925.15]
    From ClinVar: VCV000234236.13
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4097G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        4097,
        4097,
        exon=11,
        intron=None,
        var_protein="p.Gly1366Asp",
    )
    variant = VariantInfo(
        "17",
        43091032,
        43091032,
        "BRCA1",
        [VARTYPE.MISSENSE_VARIANT],
        "C",
        "T",
    )
    return (transcript, variant)


def create_example_mis_BRCA2_1() -> tuple:
    """
    Create BRCA2 missense SNV variant
    + strand, 2 position in codon
    From ClinVar: VCV000924286.2
    """
    hgvs = hgvs_parser.parse_c_posedit("c.6886A>T".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        6886,
        6886,
        exon=12,
        intron=None,
        var_protein="p.Ile2296Leu",
    )
    variant = VariantInfo(
        "13",
        32344602,
        32344602,
        "BRCA2",
        [VARTYPE.MISSENSE_VARIANT],
        "A",
        "T",
    )
    return (transcript, variant)


def create_example_mis_BRCA2_2() -> tuple:
    """
    Create BRCA2 missense SNV variant
    + strand, 2 position in codon
    From ClinVar: VCV001755970.1
    """
    hgvs = hgvs_parser.parse_c_posedit("c.6887T>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        6887,
        6887,
        exon=12,
        intron=None,
        var_protein="p.Ile2296Lys",
    )
    variant = VariantInfo(
        "13",
        32344603,
        32344603,
        "BRCA2",
        [VARTYPE.MISSENSE_VARIANT],
        "T",
        "A",
    )
    return (transcript, variant)


def create_example_mis_BRCA2_3() -> tuple:
    """
    Create BRCA2 missense SNV variant
    + strand, 2 position in codon
    From ClinVar: VCV000252458.16
    """
    hgvs = hgvs_parser.parse_c_posedit("c.6888A>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        6888,
        6888,
        exon=12,
        intron=None,
        var_protein="p.Ile2296Met",
    )
    variant = VariantInfo(
        "13",
        32344604,
        32344604,
        "BRCA2",
        [VARTYPE.MISSENSE_VARIANT],
        "A",
        "G",
    )
    return (transcript, variant)


def create_example_mis_BRCA2_2_exons() -> tuple:
    """
    Create BRCA2 missense SNV variant
    + strand, 2 position in codon
    Codon spans multiple exons
    Codon genomic pos: [VCV000038130.19, VCV000038131.38, VCV000643221.6]
    ClinVar ids: [32936829, 329368630, 32937316]
    From ClinVar: VCV000038131.38
    """
    hgvs = hgvs_parser.parse_c_posedit("c.7976G>A".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        [VARTYPE.MISSENSE_VARIANT],
        hgvs,
        7976,
        7976,
        exon=12,
        intron=None,
        var_protein="p.Arg2659Lys",
    )
    variant = VariantInfo(
        "13",
        32362693,
        32362693,
        "BRCA2",
        [VARTYPE.MISSENSE_VARIANT],
        "G",
        "A",
    )
    return (transcript, variant)
