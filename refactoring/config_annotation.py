#!/usr/bin/env python3

import yaml
import pathlib
from jsonschema import validate
from functools import reduce
import operator as op
from dataclasses import dataclass
from typing import Callable, Union, Any
from enum import Enum
from functools import partial

import refactoring.acmg_rules as rules
from refactoring.clinvar_annot import get_annotate_clinvar
from refactoring.information import (
    classification_information,
    classification_information_groups,
)
from refactoring.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
    annotate_transcripts,
)
from refactoring.var_type import VARTYPE_GROUPS
from refactoring.variant import Variant


def perform_annotation(path_config: pathlib.Path, variant: Variant):
    """
    Handle all configuration and annotation steps
    """
    config = import_config(path_config)
    final_config = select_gene_specific_functions_if_applicable(
        config, variant.variant_info.gene_name
    )
    annotations_needed = from_rule_list_get_annotations_needed(final_config["rules"])
    annotations_set_up = get_annotation_functions(
        list(annotations_needed), variant, final_config
    )
    annotations = execute_annotation(annotations_set_up)
    return annotations


def import_config(path_config: pathlib.Path) -> dict:
    """
    Import configuration and validate it
    """
    with open(path_config) as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    if not validate_yaml(config):
        raise ValueError(
            "YAML configuratino could not be validated. Plwase recheck YAML"
        )
    return config


def validate_yaml(config: dict) -> bool:
    """
    Validate yaml using a predefine json schema
    """
    json_schema_path = pathlib.Path(
        "/home/katzkean/variant_classification/refactoring/config_schema.json"
    )
    with open(json_schema_path) as f:
        json_schema = yaml.load(f, Loader=yaml.SafeLoader)
    try:
        validate(config, json_schema)
    except:
        return False
    return True


def select_gene_specific_functions_if_applicable(config: dict, gene_name: str) -> dict:
    """
    Check if gene_specific_configs are available for gene_name
    If available return gene specific configuration otherwise return standard configuration
    """
    if "gene_specific_configs" in config.keys():
        if gene_name.lower in config["gene_specific_configs"].keys():
            dir_gene_config = pathlib.Path(config["gene_specific_configs"]["root"])
            file_gene_config = config["gene_specific_configs"][gene_name.lower]
            path_gene_config = dir_gene_config / file_gene_config
            gene_config = import_config(path_gene_config)
            return gene_config
        else:
            return config
    else:
        return config


def from_rule_list_get_annotations_needed(
    rule_list: list[str],
) -> set[classification_information]:
    """
    Based on rule get classification_information objects required to apply the rules
    """
    RULE_DICTIONARY = {
        "pvs1": rules.pvs1,
        "ps1": rules.ps1,
        "pm1": rules.pm1,
        "pm2": rules.pm2,
        "pm4": rules.pm4,
        "pm5": rules.pm5,
        "pp3": rules.bp4_pp3,
        "ba1": rules.ba1,
        "bs1": rules.bs1,
        "bp3": rules.bp3,
        "bp4": rules.bp4_pp3,
        "bp7": rules.bp7,
    }
    annotations_needed = set()
    for rule in rule_list:
        try:
            rule_class = RULE_DICTIONARY[rule.lower()]
        except:
            ValueError(
                f"{rule.upper()} not valid rule. Valid rules are {RULE_DICTIONARY.keys()}"
            )
            break
        annotations_needed.update(rule_class.arguments)
    return annotations_needed



def get_annotation_functions(
    annotations_needed: list[classification_information], variant: Variant, config: dict
) -> list[classification_information]:
    """
    Based on needed annotations perform annotation
    """
    ### Dictionary for all classification_information objects that are defined in the Variant object
    ### For definition Variant object see variant.py
    dict_annotation_variant = {
        classification_information.VARIANT: partial(lambda variant: variant.variant_info, variant),
        classification_information.TRANSCRIPT: partial(lambda variant, : variant.transcript_info, variant),
        classification_information.VARIANT_HOTSPOT: partial(lambda variant : variant.affected_region.critical_region, variant),
        classification_information.VARIANT_GNOMAD: partial(lambda variant : variant.gnomad, variant),
        classification_information.VARIANT_FLOSSIES: partial(lambda variant : variant.flossies, variant),
        classification_information.VARIANT_PREDICTION: partial(lambda variant : variant.prediction_tools, variant),
    }

    ### Dictionary for all classification_information objects that have a get_annotation_function
    dict_annotation = {
        classification_information.ANNOTATED_TRANSCRIPT_LIST: lambda variant, config: get_annotation_function_annotated_transcript(
            variant, config
        ),
        classification_information.ANNOTATED_TRANSCRIPT_LIST_ACMG_Spec: lambda variant, config: get_annotation_function_annotated_transcript_acmg(
            variant, config
        ),
        classification_information.VARIANT_CLINVAR: lambda variant, config: get_annotation_function_variant_clinvar(
            variant, config
        ),
    }

    ### Dictionary for all classification_information that belong to a classification_information_groups
    dict_annotation_groups = {
        classification_information_groups.PATH: lambda annot, config: partial(
            get_path_from_config, annot.value.config_location, config
        ),
        classification_information_groups.THRESHOLDS_SINGLE: lambda annot, config: partial(
            get_threshold_from_config, annot.value.config_location, config
        ),
        classification_information_groups.THRESHOLD_PREDICTION: lambda annot, config: partial(
            get_threshold_prediction_from_config, annot.value.config_location, config
        ),
    }

    for annotation in annotations_needed:
        if annotation in dict_annotation_variant.keys():
            annotation.value.compute_function = dict_annotation_variant[annotation]
        elif annotation in dict_annotation.keys():
            annotation.value.compute_function = dict_annotation[annotation](
                variant, config
            )
        else:
            if annotation in classification_information_groups.PATH.value:
                annotation.value.compute_function = dict_annotation_groups[
                    classification_information_groups.PATH
                ](annotation, config)
            elif (
                annotation in classification_information_groups.THRESHOLDS_SINGLE.value
            ):
                annotation.value.compute_function = dict_annotation_groups[
                    classification_information_groups.THRESHOLDS_SINGLE
                ](annotation, config)
            elif (
                annotation
                in classification_information_groups.THRESHOLD_PREDICTION.value
            ):
                annotation.value.compute_function = dict_annotation_groups[
                    classification_information_groups.THRESHOLD_PREDICTION
                ](annotation, config)
            else:
                raise ValueError(f"No annotation function defined for {annotation}.")
    return annotations_needed


def get_annotation_function_annotated_transcript(
    variant: Variant, config: dict
) -> Callable[[], Any]:
    """
    Create annotation function for construction of classification_information.ANNOTATED_TRANSCIPT_LIST
    """
    relevant_classes = {
        VARTYPE_GROUPS.EXONIC: TranscriptInfo_exonic,
        VARTYPE_GROUPS.INTRONIC: TranscriptInfo_intronic,
        VARTYPE_GROUPS.START_LOST: TranscriptInfo_start_loss,
    }
    fun_dict = {}
    for name, entry in relevant_classes.items():
        fun_annot = prepare_function_for_annotation(entry.get_annotate, variant, config)
        fun_dict[name] = fun_annot
    fun = partial(annotate_transcripts, variant, fun_dict)
    return fun


def get_annotation_function_annotated_transcript_acmg(
    variant: Variant, config: dict
) -> Callable[[], Any]:
    """
    Create annotation function for construction of classification_information.ANNOTATED_TRANSCIPT_LIST_ACMG_Spec
    """
    relevant_classes = {
        VARTYPE_GROUPS.EXONIC: TranscriptInfo_exonic,
        VARTYPE_GROUPS.INTRONIC: TranscriptInfo_intronic,
        VARTYPE_GROUPS.START_LOST: TranscriptInfo_start_loss,
    }
    fun_dict = {}
    for name, entry in relevant_classes.items():
        fun_annot = prepare_function_for_annotation(entry.get_annotate, variant, config)
        fun_dict[name]["general"] = fun_annot
        fun_acmg = prepare_function_for_annotation(
            entry.get_annotation_acmg, variant, config
        )
        fun_dict[name]["acmg"] = fun_acmg
    fun = partial(annotate_transcripts, variant, fun_dict)
    return fun


def get_annotation_function_variant_clinvar(
    variant: Variant, config: dict
) -> Callable[[], Any]:
    """
    Create annotation function for construction of classification_information.VARIANT_CLINVAR
    """
    fun = prepare_function_for_annotation(get_annotate_clinvar, variant, config)
    return fun


def prepare_function_for_annotation(
    fun: Callable[[], tuple[Callable, tuple[classification_information, ...]]],
    variant: Variant,
    config: dict,
) -> Callable[[], Any]:
    """
    Prepare annotation function
    """
    annot_fun, annot_fun_args = fun()
    set_args = get_annotation_functions(list(annot_fun_args), variant, config)
    args = execute_annotation(set_args)
    fun_annot = partial(annot_fun, *[argument.value.value for argument in args])
    return fun_annot


def execute_annotation(
    annotations_to_execute: list[classification_information],
) -> list[classification_information]:
    """
    Perform annotation for entrys in list
    """
    for annotation in annotations_to_execute:
        if annotation.value.compute_function is not None:
            try:
                annotation.value.value = annotation.value.compute_function()
            except Exception:
                raise ValueError(
                    f"{annotation} function does not execute. Please check that all necessary information is available"
                )
    return annotations_to_execute


def get_path_from_config(
    config_location: Union[tuple[str, ...], None], config: dict
) -> pathlib.Path:
    """
    From config reconstruct full path
    """
    if config_location is not None:
        root_files = pathlib.Path(config[config_location[0]]["root"])
        dir_files = (
            root_files
            / pathlib.Path(config[config_location[0]][config_location[1]]["root"])
            / pathlib.Path(config[config_location[0]][config_location[1]]["version"])
        )
        file_path = dir_files / pathlib.Path(
            config[config_location[0]][config_location[1]][config_location[2]]
        )
        return file_path
    else:
        raise ValueError(
            "Accessing configuration not possible, location in configuration not specified"
        )


def get_threshold_from_config(
    config_location: Union[tuple[str, ...], None], config: dict
) -> float:
    """
    Get threshold from config
    """
    if config_location is not None:
        config_value = reduce(op.getitem, config_location, config)
        return float(config_value)
    else:
        raise ValueError(
            "Accessing configuration not possible, location in configuration not specified."
        )


def get_threshold_prediction_from_config(
    config_location: Union[tuple[str, ...], None], config: dict
):
    """
    Get threshold for prediction tools
    """
    if config_location is not None:
        config_prediction_tool = reduce(op.getitem, config_location, config)
        name = config_prediction_tool["name"]
        thr_benign = config_prediction_tool["benign"]
        thr_pathogenic = config_prediction_tool["pathogenic"]
        if thr_benign < thr_pathogenic:
            return Two_threshold(
                name,
                thr_benign,
                THRESHOLD_DIRECTION.LOWER,
                thr_pathogenic,
                THRESHOLD_DIRECTION.HIGHER,
            )
        else:
            return Two_threshold(
                name,
                thr_benign,
                THRESHOLD_DIRECTION.HIGHER,
                thr_pathogenic,
                THRESHOLD_DIRECTION.LOWER,
            )
    else:
        raise ValueError(
            "Accessing configuration not possible, location in configuration not specified."
        )


@dataclass
class Threshold:
    name: str


class THRESHOLD_DIRECTION(Enum):
    HIGHER = "higher"
    LOWER = "lower"


@dataclass
class One_threshold(Threshold):
    threshold: float
    direction: THRESHOLD_DIRECTION


@dataclass
class Two_threshold(Threshold):
    threshold_benign: float
    direction_benign: THRESHOLD_DIRECTION
    threshold_pathogenic: float
    direction_pathogenic: THRESHOLD_DIRECTION
