#!/usr/bin/env python3

from variant import VariantInfo, TranscriptInfo

import hgvs.parser
import hgvs.posedit

hgvs_parser = hgvs.parser.Parser()


def create_example_splice_acceptor() -> tuple(TranscriptInfo, VariantInfo):
    """
    Create BR
    """
    hgvs = hgvs_parser.parse_c_posedit("c.5278-1G>A".split("c.")[1])
    transcript = TranscriptInfo("ENST00000357654", "splice_acceptor", hgvs, 5278, 5278)
    variant = VariantInfo(
        "BRCA1", ["splice_acceptor"], "17", 41203135, 41203135, "some_id", "A", "G"
    )
    return (transcript, variant)


def create_example_splice_donor() -> tuple(TranscriptInfo, VariantInfo):
    """
    Create BR
    """
    hgvs = hgvs_parser.parse_c_posedit("c.4986+1G>T".split("c.")[1])
    transcript = TranscriptInfo("ENST00000357654", "splice_donor", hgvs, 4986, 4986)
    variant = VariantInfo(
        "BRCA1", ["splice_donor"], "17", 41222944, 41222944, "some_id", "G", "A"
    )
    return (transcript, variant)
