#!/usr/bin/env python3

from variant_classification.variant import (
    PopulationDatabases_gnomAD,
    Variant,
    PopulationDatabases,
    AffectedRegion,
)
from test.example_variant_splicing_GRCh38 import create_example_splice_acceptor_BRCA1

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
