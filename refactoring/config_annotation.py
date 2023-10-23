#!/usr/bin/env python3

from cyvcf2 import Variant
import yaml
import pathlib
from jsonschema import validate
from functools import reduce
import operator as op
from dataclasses import dataclass
from enum import Enum
from functools import partial

import refactoring.acmg_rules as rules
from refactoring.information import classification_information, classification_information_groups
from refactoring.variant import Variant


def perform_annotation(path_config: pathlib.Path, variant: Variant):
    config = import_config(path_config)
    final_config = select_gene_specific_functions_if_applicable(config, variant.variant_info.gene_name)
    annotations_needed = from_rule_list_get_annotations_needed(final_config["rules"])
    return annotations_needed


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
        print(rule)
        try:
            rule_class = RULE_DICTIONARY[rule.lower()]
        except:
            print(
                f"{rule.upper()} not valid rule. Valid rules are {RULE_DICTIONARY.keys()}"
            )
            break
        annotations_needed.update(rule_class.arguments)
    return annotations_needed



def set_up_annotations(annotations_needed: list[classification_information], variant: Variant, config: dict):
    """
    Based on needed annotations perform annotation
    """
    for annotation in annotations_needed:
        # Update python to current version in order for match to work
        if annotation is classification_information.VARIANT:
            classification_information.VARIANT.value.value = variant.variant_info
        elif annotation is classification_information.TRANSCRIPT:
            print("Hello from Transcript")
            classification_information.TRANSCRIPT.value.value = variant.transcript_info
        elif annotation is classification_information.ANNOTATED_TRANSCRIPT_LIST:
            print("Hello from Annotated transcript")
            fun = lambda: ""
        elif annotation is classification_information.ANNOTATED_TRANSCRIPT_LIST_ACMG_Spec:
            print("Hello from ACMG transcript")
            fun = lambda: ""
        elif classification_information.VARIANT_CLINVAR:
            print("Hello from ClinVar")
            fun = lambda: ""
        elif annotation in classification_information_groups.PATH.value:
            print("Hello from Path")
            fun = partial(create_path_from_config, annotation.value.config_location, config)
        elif annotation in classification_information_groups.THRESHOLDS_SINGLE.value:
            print("This is threshold BA1")
            fun = partial(get_threshold_from_config, annotation.value.config_location, config)
        elif annotation in classification_information_groups.THRESHOLD_PREDICTION.value:
            fun = parital()


def create_path_from_config(config_location: tuple[str, ...], config: dict) -> pathlib.Path:
    """
    From config reconstruct full path
    """
    if config_location is not None:
        root_files = pathlib.Path(config[config_location[0]]["root"])
        dir_files = root_files / pathlib.Path(config[config_location[0]][config_location[1]]) / pathlib.Path(config[config_location[0]][config[config_location[1]]]["version"])
        file_path = dir_files / pathlib.Path(config[config_location[0]][config_location[1]][config_location[2]])
        return file_path
    else:
        raise ValueError("Accessing configuration not possible, location in configuration not specified")


def get_threshold_from_config(config_location: tuple[str, ...], config: dict) -> float:
    """
    Get threshold from config
    """
    if config_location is not None:
        config_value = reduce(op.getitem, config_location , config)
        return float(config_value)
    else:
        raise ValueError("Accessing configuration not possible, location in configuration not specified.")


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
