#!/usr/bin/env python3

import yaml
import pathlib
from jsonschema import validate


def load_config(path_config: pathlib.Path) -> dict:
    """
    Import configuration and validate it
    """
    with open(path_config) as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    if not validate_config(config):
        raise ValueError(
            "YAML configuration could not be validated. Please recheck YAML"
        )
    return config


def validate_config(config: dict) -> bool:
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
    except Exception:
        return False
    return True


def get_gene_specific_config(config: dict, gene_name: str) -> dict:
    """
    Check if gene_specific_configs are available for gene_name
    If available return gene specific configuration otherwise return standard configuration
    """
    if "gene_specific_configs" in config.keys():
        if gene_name.lower in config["gene_specific_configs"].keys():
            dir_gene_config = pathlib.Path(config["gene_specific_configs"]["root"])
            file_gene_config = config["gene_specific_configs"][gene_name.lower]
            path_gene_config = dir_gene_config / file_gene_config
            gene_config = load_config(path_gene_config)
            return gene_config
        else:
            return config
    else:
        return config
