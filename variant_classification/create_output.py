#!/usr/bin/env python3

import json
import pathlib

from jsonschema import validate

from variant_classification.acmg_rules.utils import RuleResult


def create_output(rule_results: list[RuleResult]) -> str:
    """
    From list of RuleResult object that meets the classified schema
    """
    out_dict = {}
    for result in rule_results:
        result_dict = result.create_dict()
        out_dict = out_dict | result_dict
    if not validate_output(out_dict):
        raise ValueError("Output could not be validated. Please check.")
    result_json = json.dumps(out_dict)
    return result_json


def validate_output(out_dict: dict) -> bool:
    """
    Validate output
    """
    json_schema_path = pathlib.Path("./API/schema_output_acmg.json")
    with open(json_schema_path) as f:
        json_schema = json.load(f)
    try:
        validate(out_dict, json_schema)
    except Exception:
        return False
    return True
