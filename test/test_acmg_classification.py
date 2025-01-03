#!/usr/bin/env python

import pandas as pd

from variant_classification.classification_schemata.utils import (
    get_final_classification_from_possible_classes,
)
from variant_classification.final_classification import get_classification

# def test_select_correct_schema():


def test_bp1_BRCA():
    test_results = {
        "BP1": {
            "rule_type": "general",
            "evidence_type": "pathogenic",
            "status": True,
            "strength": "strong",
            "comment": "test",
        },
        "PS1": {
            "rule_type": "protein",
            "evidence_type": "pathogenic",
            "status": False,
            "strength": "strong",
            "comment": "Test",
        },
    }
    rules_df = pd.DataFrame(test_results).transpose()
    applicable_rules = rules_df[rules_df.status == True]
    config_name = "ACMG BRCA1"
    version = "1.1.0"
    class_result = get_classification(applicable_rules, config_name, version)
    assert class_result == 2


def test_bp1_non_BRCA():
    test_results = {
        "BP1": {
            "rule_type": "general",
            "evidence_type": "pathogenic",
            "status": True,
            "strength": "strong",
            "comment": "test",
        },
        "PS1": {
            "rule_type": "protein",
            "evidence_type": "pathogenic",
            "status": False,
            "strength": "strong",
            "comment": "Test",
        },
    }
    rules_df = pd.DataFrame(test_results).transpose()
    applicable_rules = rules_df[rules_df.status == True]
    config_name = "ACMG PALB2"
    version = "1.1.0"
    class_result = get_classification(applicable_rules, config_name, version)
    assert class_result == 3


def test_lb_and_b_to_b():
    pos_class = [1, 2]
    final_class = get_final_classification_from_possible_classes(pos_class)
    assert final_class == 1


def test_lp_and_p_to_p():
    pos_class = [4, 5]
    final_class = get_final_classification_from_possible_classes(pos_class)
    assert final_class == 5


def test_lp_and_lb_to_vus():
    pos_class = [2, 4]
    final_class = get_final_classification_from_possible_classes(pos_class)
    assert final_class == 3


def test_p_and_b_to_vus():
    pos_class = [1, 5]
    final_class = get_final_classification_from_possible_classes(pos_class)
    assert final_class == 3


def test_p_and_p_to_p():
    pos_class = [5, 5]
    final_class = get_final_classification_from_possible_classes(pos_class)
    assert final_class == 5


def test_b_and_b_to_b():
    pos_class = [1, 1]
    final_class = get_final_classification_from_possible_classes(pos_class)
    assert final_class == 1


def test_no_class_to_vus():
    pos_class = []
    final_class = get_final_classification_from_possible_classes(pos_class)
    assert final_class == 3
