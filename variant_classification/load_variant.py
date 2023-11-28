#!/usr/bin/env python3

import pathlib
import json
import logging

from typing import Optional
from jsonschema import validate
import hgvs.parser
import hgvs.posedit
from var_type import VARTYPE
from variant import (
    FunctionalData,
    MultifactorialLikelihood,
    Variant,
    VariantInfo,
    TranscriptInfo,
    PopulationDatabases_gnomAD,
    PopulationDatabases,
    AffectedRegion,
)
from os import path


logger = logging.getLogger("HerediClass.load_variant")


hgvs_parser = hgvs.parser.Parser()


def load_variant(var_str: str) -> Variant:
    """
    Load variant from json string
    """
    var_dict = json.loads(var_str)
    if not validate_variant(var_dict):
        raise ValueError("Variant could not be validated. Please check.")
    variant = create_variant(var_dict)
    return variant


def validate_variant(var_dict: dict) -> bool:
    """
    Validate variant input
    """
    base_path =  path.dirname(path.dirname(path.abspath(__file__)))
    json_schema_path = pathlib.Path(path.join(base_path, "API/schema_input.json"))
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
        affected_region=affected_region,
    )
    return variant


def create_variantInfo(variant_json: dict) -> VariantInfo:
    """
    Create VariantInfo object from variant_json
    """
    chr = variant_json["chr"]
    var_type = get_vartype_list(variant_json["variant_type"])
    gene_name = variant_json["gene"]
    ref = variant_json["ref"]
    alt = variant_json["alt"]
    if len(ref) == 1 and len(alt) == 1:
        # Create genomic positon for substitution
        genomic_start = variant_json["pos"]
        genomic_end = variant_json["pos"]
    elif len(ref) > len(alt):
        # Create genomic position for deletions
        genomic_start = variant_json["pos"] + 1
        del_length = len(alt) - 1
        genomic_end = genomic_start + del_length - 1
        alt = ""
        ref = ref[1:]
    elif len(ref) < len(alt):
        # Create genomic position for insertions
        genomic_start = variant_json["pos"]
        genomic_end = variant_json["pos"]
        ref = ""
        alt = alt[1:]
    else:
        raise ValueError(
            f"Variant is not of type substitution, insertion or deletion. Please check variant input."
        )
    var_info = VariantInfo(
        chr=chr,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
        gene_name=gene_name,
        var_type=var_type,
        var_ref=ref,
        var_obs=alt,
    )
    return var_info


def create_transcriptinfo(variant_json: dict) -> list[TranscriptInfo]:
    """
    Create TranscriptInfo object from variant_json
    """
    transcripts_dict = variant_json["variant_effect"]
    transcripts = []
    for trans_dict in transcripts_dict:
        transcript_id = trans_dict["transcript"]
        hgvs_c_str = trans_dict["hgvs_c"]
        hgvs_c = hgvs_parser.parse_c_posedit(hgvs_c_str.split("c.")[1])
        var_start = hgvs_c.pos.start.base
        var_stop = hgvs_c.pos.end.base
        var_type = get_vartype_list(variant_json["variant_type"])
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
            var_stop=var_stop,
            var_protein=hgvs_p,
            exon=exon,
            intron=intron,
        )
        transcripts.append(transcript)
    return transcripts


def get_vartype_list(var_type_str: list[str]) -> list[VARTYPE]:
    """
    From list of var_types in str format produce list of VARTYPE Enums
    """
    var_types = []
    for entry in var_type_str:
        try:
            var_type = VARTYPE(entry)
            var_types.append(var_type)
        except ValueError:
            continue
    if len(var_types) == 0:
        raise ValueError(
            f"For the variant types {var_type_str} no entry in VARTYPE could be found. Please check."
        )
    return var_types


def create_prediction_dict(variant_json: dict) -> Optional[dict[str, float]]:
    """
    Create predciton tool dictionary from variant_json
    """
    try:
        patho_prediction = variant_json["pathogenicity_prediction_tools"]
    except KeyError:
        patho_prediction = {}
    try:
        splice_prediction = variant_json["splicing_prediction_tools"]
    except KeyError:
        splice_prediction = {}
    prediction = patho_prediction | splice_prediction
    if not bool(prediction):
        return None
    return prediction


def create_gnomad(variant_json: dict) -> Optional[PopulationDatabases_gnomAD]:
    """
    Create gnomAD object from variant_json
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
        name=name,
        frequency=frequency,
        count=allele_count,
        popmax=popmax,
        popmax_frequency=popmax_AF,
        popmax_allele_count=popmax_AC,
    )
    return gnomad


def create_flossies(variant_json: dict) -> Optional[PopulationDatabases]:
    """
    Create FLOSSIES object from variant_json
    """
    try:
        flossies_dict = variant_json["FLOSSIES"]
    except KeyError:
        return None
    name = "flossies"
    if flossies_dict["AFR"] > flossies_dict["EUR"]:
        count = flossies_dict["AFR"]
    else:
        count = flossies_dict["EUR"]
    return PopulationDatabases(name=name, count=count, frequency=None)


def create_affected_region(variant_json: dict) -> Optional[AffectedRegion]:
    """
    Create AffectedRegion object from variant_json
    """
    try:
        crit_region = variant_json["VUS_task_force_domain"]
    except KeyError:
        crit_region = False
    try:
        cancer_hotspot = variant_json["cancer_hotspot"]
    except KeyError:
        cancer_hotspot = False
    try:
        cold_spot = variant_json["cold_spot"]
    except KeyError:
        cold_spot = False
    aff_reg = AffectedRegion(
        critical_region=crit_region, cancer_hotspot=cancer_hotspot, cold_spot=cold_spot
    )
    return aff_reg


def get_mutlifactorial_likelihood(
    variant_json: dict,
) -> Optional[MultifactorialLikelihood]:
    """
    Create MultifactorialLikelihood object from variant_json
    """
    try:
        prior = variant_json["prior"]
    except KeyError:
        prior = None
    try:
        co_occurrence = variant_json["co-occurrence"]
    except KeyError:
        co_occurrence = None
    try:
        segregation = variant_json["segregation"]
    except KeyError:
        segregation = None
    try:
        multifactorial_likelihood = variant_json["multifactorial_log-likelihood"]
    except KeyError:
        multifactorial_likelihood = None
    if not any([prior, co_occurrence, segregation, multifactorial_likelihood]):
        return None
    multfaclike = MultifactorialLikelihood(
        prior=prior,
        multifactorial_likelihood=multifactorial_likelihood,
        co_occurrence=co_occurrence,
        co_segregation=segregation,
    )
    return multfaclike


def create_functional_data(key: str, variant_json: dict) -> Optional[FunctionalData]:
    """
    Create FunctionalData object from variant_json
    """
    try:
        func_data = variant_json[key]
    except KeyError:
        return None
    patho = func_data["pathogenic"]
    ben = func_data["benign"]
    performed = func_data["performed"]
    func = FunctionalData(performed=performed, pathogenic=patho, benign=ben)
    return func
