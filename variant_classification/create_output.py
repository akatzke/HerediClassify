#!/usr/bin/env python3

import json

from variant_classification.acmg_rules.utils import RuleResult


def create_output(rule_results: list[RuleResult]) -> str:
    """
    From list of RuleResult object that meets the classified schema
    """
    out_dict = {}
    for result in rule_results:
        result_dict = result.create_dict()
        out_dict = out_dict | result_dict
    result_json = json.dumps(out_dict)
    return result_json
