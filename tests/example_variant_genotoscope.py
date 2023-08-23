#!/usr/bin/env python3

from refactoring.variant import VariantInfo, TranscriptInfo

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_brca2_missense() -> tuple:
    """
    Variant defined in origianl GenOtoScope test set
    """
    hgvs = hgvs_parser.parse_c_posedit("c.3073A>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000380152",
        ["missense_variant"],
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
        ["frameshift_variant"],
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
        ["frameshift_variant"],
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
