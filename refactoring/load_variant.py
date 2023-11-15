#!/usr/bin/env python3

import pathlib
import json
import io

from typing import Optional
from jsonschema import validate
from refactoring.load_config import validate_config

from refactoring.variant import Variant, VariantInfo, TranscriptInfo, PopulationDatabases_gnomAD, PopulationDatabases, AffectedRegion

def load_variant(var_str: str) -> Variant:
    """
    Load variant from json string
    """
    var_dict = json.load(var_str)
    if not validate_variant(var_dict):
        raise ValueError(
            "Variant could not be validated. Please check."
            )
    variant = create_variant(var_dict)
    return variant


def validate_variant(var_dict: dict) -> bool:
    """
    Validate variant input
    """
    json_schema_path = pathlib.Path("/home/katzkean/variant_classification/API/schema_input.json")
    with open(json_schema_path) as f:
        json_schema = json.load(f)
    try:
        validate(var_dict, json_schema)
    except Exception:
        return False
    return True


def create_variant(variant_json: dict) -> Variant:
    """
    Create Variant object from variant_json
    """
    var_info = create_variantInfo(variant_json)
    trans_info_list = create_transcriptinfo(variant_json)
    prediction_tools = create_prediction_dict(variant_json)
    gnomad = create_gnomad(variant_json)
    flossies = create_flossies(variant_json)
    affected_region = create_affected_region(variant_json)
    variant = Variant(
        variant_info=var_info,
        transcript_info=trans_info_list,
        prediction_tools=prediction_tools,
        gnomad=gnomad,
        flossies=flossies,
        affected_region=affected_region
    )
    return variant


def create_variantInfo(variant_json: dict) -> VariantInfo:
    """
    Create VariantInfo object from variant_json
    """
    chr = variant_json["chr"]
    genomic_start = variant_json["pos"]
    genomic_end = variant_json["pos"]
    var_type =  variant_json["variant_type"]
    gene_name = variant_json["gene"]
    ref = variant_json["ref"]
    alt = variant_json["alt"]
    var_info = VariantInfo(
        chr = chr,
        genomic_start = genomic_start,
        genomic_end=genomic_end,
        gene_name=gene_name,
        var_type=var_type,
        var_ref = ref,
        var_obs = alt,
        )
    return var_info


def create_transcriptinfo(variant_json: dict) -> list[TranscriptInfo]:
    """
    Create TranscriptInfo object from variant_json
    """
    transcripts_dict = variant_json["variant_effect"]:
    transcripts = []
    for trans_dict in transcripts_dict:
        transcript_id = trans_dict["transcript"]
        hgvs_c_str = trans_dict["hgvs_c"]
        hgvs_c = hgvs_parser.parse_c_posedit(hgvs_c_str.split("c.")[1])
        var_start = hgvs_c.pos.start.base
        var_stop = hgvs_c.pos.end.base
        var_type = trans_dict["variant_type"]
        try:
            hgvs_p = trans_dict["hgvs_p"]
        except KeyError:
            hgvs_p = None
        try:
            exon = trans_dict["exon"]
        except KeyError:
            exon = None
        try:
            intron = trans_dict["intron"]
        except:
            intron = None
        transcript = TranscriptInfo(
            transcript_id=transcript_id,
            var_type=var_type,
            var_hgvs=hgvs_c,
            var_start=var_start,
            var_stop = var_stop,
            var_protein=hgvs_p,
            exon = exon,
            intron = intron
            )
        transcripts.append(transcript)
    return transcripts

def create_prediction_dict(variant_json: dict) -> Optional[dict[str, float]]:
    """
    Create predciton tool dictionary
    """
    patho_prediction = variant_json["pathogenicity_prediction_tools"]
    splice_prediction = variant_json["splicing_prediction_tools"]
    prediction = patho_prediction | splice_prediction
    return prediction

def create_gnomad(variant_json: dict) -> Optional[PopulationDatabases_gnomAD]:
    """
    Create gnomAD object
    """
    try:
        gnomad_dict = variant_json["gnomAD"]
    except KeyError:
        return None
    name = "gnomAD"
    frequency = gnomad_dict["AF"]
    allele_count = gnomad_dict["AC"]
    popmax = gnomad_dict["popmax"]
    popmax_AF = gnomad_dict["popmax_AF"]
    popmax_AC = gnomad_dict["popmax_AC"]
    gnomad = PopulationDatabases_gnomAD(
        name = name,
        frequency=frequency,
        allele_count=allele_count,
        popmax=popmax,
        popmax_frequency=popmax_AF,
        popmax_allele_count=popmax_AC
    )
    return gnomad

def create_flossies(variant_json: dict) -> Optional[PopulationDatabases]:
    try:
        flossies_dict = variant_json["FLOSSIES"]
    except KeyError:
        return None
    name = "flossies"
    if flossies_dict["AFR"] > flossies_dict["EUR"]:
        count = flossies_dict["AFR"]
    else:
        count = flossies_dict["EUR"]
    return PopulationDatabases(
        name = name,
        frequency=count
    )

def create_affected_region(variant_json: dict) -> Optional[AffectedRegion]:
    crit_region = variant_json["VUS_task_force_domain"]
    cancer_hotspot = variant_json["cancer_hotspot"]
    cold_spot = variant_json["cold_spot"]
    aff_reg = AffectedRegion(
        critical_region= crit_region,
        cancer_hotspot=cancer_hotspot,
        cold_spot=cold_spot
        )
    return aff_reg
