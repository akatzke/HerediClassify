#!/usr/bin/env python3

from collections.abc import Callable, Collection
import yaml
import pathlib
from dataclasses import dataclass
from jsonschema import validate
from enum import Enum

import refactoring.acmg_rules as rules
import refactoring.information as info
from refactoring.transcript_annotated import (
    annotate_transcripts,
    annotate_transcripts_acmg_specification,
)
from refactoring.clinvar_annot import annotate_clinvar
from refactoring.variant import VariantInfo


def import_config(path_config: pathlib.Path) -> dict:
    """
    Import configuration
    """
    with open(path_config) as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    if not validate_yaml(config):
        raise ValueError(
            "YAML configuratino could not be validated. Plwase recheck YAML"
        )
    return config


def create_complete_file_paths(config: dict) -> dict:
    """
    Paths are given in roots and subfolders etc
    Here should be documented how to create complete paths from these partial files
    """
    return config


def select_gene_specific_functions_if_applicable(config: dict, variant: VariantInfo):
    """
    Check if config has gene_specific_annotations_defined
    If not defined simple use standrad config
    If gene specific configs are listed here, check if gene has gene_specific config
    Replace config with gene_specific config
    """
    pass


def get_annotatins(config: dict):
    """
    Based on the annotation and loaded variant information
    This functions
    """
    annotations_needed = from_rule_list_get_annotations_needed(config["rules"])
    for annotation in annotations_needed:
        # Update python to current version in order for match to work
        fun = None
        match annotation:
            case info.classification_information.ANNOTATED_TRANSCRIPT_LIST:
                ...
                # here goes partial definition of calulating function
                fun = lambda: ""

        annotation.value.compute_function = fun


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


def from_rule_list_get_annotations_needed(
    rule_list: list[str],
) -> set[info.classification_information]:
    RULE_DICTIONARY: dict[str, rules.abstract_rule] = {
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
            print(
                f"{rule.upper()} not valid rule. Valid rules are {RULE_DICTIONARY.keys()}"
            )
            break
        annotations_needed.update(rule_class.arguments)
    return annotations_needed


def from_annotations_list_get_annotation_function(
    annotation_list: Collection[info.classification_information],
) -> set[str]:
    ANNOTATION_DICTIONARY: dict[info.classification_information, Callable] = {
        info.classification_information.ANNOTATED_TRANSCRIPT_LIST: annotate_transcripts,
        info.classification_information.ANNOTATED_TRANSCRIPT_LIST_ACMG_Spec: annotate_transcripts_acmg_specification,
        info.classification_information.VARIANT_CLINVAR: annotate_clinvar,
    }
    annot_information_list = set()
    for annot in annotation_list:
        try:
            annot_info = ANNOTATION_DICTIONARY[annot]
        except:
            continue
        annot_information_list.add(annot_info)
    return annot_information_list


@dataclass
class Threshold:
    name: str


class THRESHOLD_DIRECTION(Enum):
    HIGHER = "higher"
    LOWER = "lower"


class One_threshold(Threshold):
    threshold: float
    direction: THRESHOLD_DIRECTION


class Two_threshold(Threshold):
    threshold_benign: float
    direction_benign: THRESHOLD_DIRECTION
    threshold_pathogenic: float
    direction_pathogenic: THRESHOLD_DIRECTION
