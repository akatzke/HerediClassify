#!/usr/bin/env python3


from variant_classification.variant import VariantInfo, TranscriptInfo, VARTYPE

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_example_snv_syn() -> tuple:
    """
    Create BRCA1 synonymous SNV variant
    Classified as class 2 (BP1, BP4, BP6, PM2_pp)
    """
    hgvs = hgvs_parser.parse_c_posedit("c.1740C>T".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000357654",
        [VARTYPE.SYNONYMOUS_VARIANT],
        hgvs,
        4986,
        4986,
        exon=10,
        intron=None,
        var_protein=None,
    )
    variant = VariantInfo(
        "17",
        43093791,
        43093791,
        "BRCA1",
        [VARTYPE.SYNONYMOUS_VARIANT],
        "G",
        "A",
    )
    return (transcript, variant)


def create_start_lost() -> tuple:
    """
    Create RAD51D start lost variant
    From HerediVar
    Classified as class 4 (PVS1, PS1, PM2)
    """
    hgvs = hgvs_parser.parse_c_posedit("c.2T>G".split("c.")[1])
    transcript = TranscriptInfo(
        "ENST00000345365",
        [VARTYPE.START_LOST],
        hgvs,
        2,
        2,
        exon=1,
        intron=None,
        var_protein="p.Met1?",
    )
    variant = VariantInfo(
        "17", 35119612, 35119612, "RAD51D", [VARTYPE.START_LOST], "A", "C"
    )
    return (transcript, variant)
