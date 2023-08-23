#!/usr/bin/env python3

from refactoring.variant import VariantInfo, TranscriptInfo

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_example_snv_syn() -> tuple:
    """
    Create BRCA1 synonymous SNV variant
    From HerediVar
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1740C>T".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        ["synonymous_variant"],
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
        ["start_lost"],
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
