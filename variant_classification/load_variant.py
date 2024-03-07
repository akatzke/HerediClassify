#!/usr/bin/env python3

import pathlib
import json
import logging

from typing import Optional
from jsonschema import validate
import hgvs.parser
import hgvs.posedit
import hgvs.exceptions
from var_type import VARTYPE
from variant import (
    ALLELIC,
    FunctionalData,
    MultifactorialLikelihood,
    RNAData,
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
    base_path = path.dirname(path.dirname(path.abspath(__file__)))
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
    gnomad_popmax = create_gnomad(variant_json, "popmax")
    gnomad_faf = create_gnomad(variant_json, "faf_popmax")
    flossies = create_flossies(variant_json)
    cancer_hotspots = create_cancer_hotspots(variant_json)
    affected_region = create_affected_region(variant_json)
    mRNA_result = create_rna_data("mRNA_analysis", variant_json)
    functional_data = create_functional_data("functional_data", variant_json)
    variant = Variant(
        variant_info=var_info,
        transcript_info=trans_info_list,
        prediction_tools=prediction_tools,
        gnomad_popmax=gnomad_popmax,
        gnomad_faf=gnomad_faf,
        flossies=flossies,
        cancerhotspots=cancer_hotspots,
        affected_region=affected_region,
        functional_assay=functional_data,
        splicing_assay=mRNA_result,
    )
    return variant


def create_variantInfo(variant_json: dict) -> VariantInfo:
    """
    Create VariantInfo object from variant_json
    """
    chr = variant_json["chr"]
    if "chr" in chr:
        chr = chr.split("chr")[1]
    var_type = get_vartype_list(variant_json["variant_type"])
    gene_name = variant_json["gene"]
    ref = variant_json["ref"]
    alt = variant_json["alt"]
    if len(ref) == 1 and len(alt) == 1:
        # Create genomic positon for substitution
        genomic_start = variant_json["pos"]
        genomic_end = variant_json["pos"]
    elif len(ref) > 1 and len(alt) > 1:
        # Create genomic position for indels
        genomic_start = variant_json["pos"] + 1
        del_length = len(ref) - 1
        genomic_end = genomic_start + del_length - 1
    elif len(ref) > 1 and len(alt) == 1:
        # Create genomic position for deletions
        genomic_start = variant_json["pos"] + 1
        del_length = len(ref) - 1
        genomic_end = genomic_start + del_length - 1
        alt = ""
        ref = ref[1:]
    elif len(ref) == 1 and len(alt) > 1:
        # Create genomic position for insertions
        genomic_start = variant_json["pos"]
        genomic_end = variant_json["pos"]
        ref = ""
        alt = alt[1:]
    else:
        # All other cases
        genomic_start = variant_json["pos"]
        genomic_end = variant_json["pos"]
    if genomic_start > genomic_end:
        raise ValueError(f"genomic_start {genomic_start} is bigger than {genomic_end}")
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
        if hgvs_c_str is None or hgvs_c_str == "None":
            continue
        if hgvs_c_str[0:2] == "n.":
            continue
        if "c.*" in hgvs_c_str:
            continue
        try:
            hgvs_c = hgvs_parser.parse_c_posedit(hgvs_c_str.split("c.")[1])
        except hgvs.exceptions.HGVSParseError:
            continue
        var_start = hgvs_c.pos.start.base
        var_stop = hgvs_c.pos.end.base
        var_type = get_vartype_list(trans_dict["variant_type"])
        hgvs_p = trans_dict.get("hgvs_p", None)
        exon = trans_dict.get("exon", None)
        intron = trans_dict.get("intron", None)
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
        entry = entry.lower().strip().replace(" ", "_")
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


def create_prediction_dict(variant_json: dict) -> dict[str, float]:
    """
    Create predciton tool dictionary from variant_json
    """
    patho_prediction = variant_json.get("pathogenicity_prediction_tools", dict())
    splice_prediction = variant_json.get("splicing_prediction_tools", dict())
    prediction = patho_prediction | splice_prediction
    return prediction


def create_gnomad(variant_json: dict, type: str) -> PopulationDatabases_gnomAD:
    """
    Create gnomAD object from variant_json
    """
    gnomad_dict = variant_json.get("gnomAD", dict())
    name = "gnomAD"
    frequency = gnomad_dict.get("AF", 0)
    allele_count = gnomad_dict.get("AC", 0)
    subpopulation = gnomad_dict.get("subpopulation", "None")
    subpopulation_AF = gnomad_dict.get(f"{type}_AF", 0)
    subpopulation_AC = gnomad_dict.get(f"popmax_AC", 0)
    gnomad = PopulationDatabases_gnomAD(
        name=name,
        frequency=frequency,
        count=allele_count,
        subpopulation=subpopulation,
        subpopulation_frequency=subpopulation_AF,
        subpopulation_allele_count=subpopulation_AC,
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


def create_cancer_hotspots(variant_json: dict) -> Optional[PopulationDatabases]:
    """
    Create cancer hotspots object from variant_json
    """
    try:
        cancer_hotspots_dict = variant_json["cancer_hotspots"]
    except KeyError:
        return None
    count = cancer_hotspots_dict.get("AC", 0)
    frequency = cancer_hotspots_dict.get("AF", 0)
    return PopulationDatabases(
        name="Cancer Hotspots",
        count=count,
        frequency=frequency,
    )


def create_affected_region(variant_json: dict) -> AffectedRegion:
    """
    Create AffectedRegion object from variant_json
    """
    crit_region = variant_json.get("VUS_task_force_domain", False)
    cold_spot = variant_json.get("cold_spot", False)
    aff_reg = AffectedRegion(critical_region=crit_region, cold_spot=cold_spot)
    return aff_reg


def get_mutlifactorial_likelihood(
    variant_json: dict,
) -> Optional[MultifactorialLikelihood]:
    """
    Create MultifactorialLikelihood object from variant_json
    """
    prior = variant_json.get("prior", None)
    co_occurrence = variant_json.get("co-occurrence", None)
    segregation = variant_json.get("segregation", None)
    multifactorial_likelihood = variant_json.get("multifactorial_log-likelihood", None)
    if not any([prior, co_occurrence, segregation, multifactorial_likelihood]):
        return None
    multfaclike = MultifactorialLikelihood(
        prior=prior,
        multifactorial_likelihood=multifactorial_likelihood,
        co_occurrence=co_occurrence,
        co_segregation=segregation,
    )
    return multfaclike


def create_functional_data(
    key: str, variant_json: dict
) -> Optional[list[FunctionalData]]:
    """
    Create FunctionalData object from variant_json
    """
    try:
        func_data = variant_json[key]
    except KeyError:
        return None
    func_list = []
    for entry in func_data:
        patho = entry["pathogenic"]
        ben = entry["benign"]
        func = FunctionalData(pathogenic=patho, benign=ben)
        func_list.append(func)
    return func_list


def create_rna_data(key: str, variant_json: dict) -> Optional[list[RNAData]]:
    """
    Create RNAData object from variant_json
    """
    try:
        func_data = variant_json[key]
    except KeyError:
        return None
    func_list = []
    for entry in func_data:
        minigene = entry["minigene"]
        patient_rna = entry["patient_rna"]
        allelic = ALLELIC(entry["allelic"].lower().strip())
        quantification = entry["quantification"]
        func = RNAData(
            minigene=minigene,
            patient_rna=patient_rna,
            allelic=allelic,
            quantification=quantification,
        )
        func_list.append(func)
    return func_list
