#!/usr/bin/env python3

import json

from test.test_import_variant import create_json_string_from_variant
import test.paths as paths

from variant_classification.classify import classify


def test_acmg_bard1_frameshift():
    """
    Test general ACMG guidelines
    """
    path_variant = paths.TEST / "test_variants" / "BARD1_frameshift_variant.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    _, results = classify(path_config, var_str)
    results_dict = json.loads(results)
    key_list = [key for key in results_dict.keys()]
    rules_apply = ["PVS1_protein", "BP4_splicing"]
    for rule in rules_apply:
        assert results_dict[rule]["status"]
    rules_not_apply = [rule for rule in key_list if rule not in rules_apply]
    for rule in rules_not_apply:
        assert not results_dict[rule]["status"]


def test_acmg_MSH6_missense_variant():
    """
    Test general ACMG guidelines
    Class 2 due to BS1_supporting
    """
    path_variant = paths.TEST / "test_variants" / "MSH6_missense_variant.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    _, results = classify(path_config, var_str)
    results_dict = json.loads(results)
    key_list = [key for key in results_dict.keys()]
    rules_apply = ["BP4_protein", "BP4_splicing"]
    for rule in rules_apply:
        assert results_dict[rule]["status"]
    rules_not_apply = [rule for rule in key_list if rule not in rules_apply]
    for rule in rules_not_apply:
        assert not results_dict[rule]["status"]


def test_acmg_MSH6_missense_variant_2():
    """
    Test general ACMG guidelines
    Class 4 due to functional study
    """
    path_variant = paths.TEST / "test_variants" / "MSH6_missense_variant_2.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    _, results = classify(path_config, var_str)
    results_dict = json.loads(results)
    key_list = [key for key in results_dict.keys()]
    rules_apply = ["PM2", "PP3_protein", "BP4_splicing"]
    for rule in rules_apply:
        assert results_dict[rule]["status"]
    rules_not_apply = [rule for rule in key_list if rule not in rules_apply]
    for rule in rules_not_apply:
        assert not results_dict[rule]["status"]


def test_acmg_PMS2_splice_donor():
    """
    Test general ACMG guidelines
    """
    path_variant = paths.TEST / "test_variants" / "PMS2_splice_donor_variant.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    _, results = classify(path_config, var_str)
    results_dict = json.loads(results)
    key_list = [key for key in results_dict.keys()]
    rules_apply = ["PVS1_splicing", "PP3_splicing", "PS1_splicing"]
    for rule in rules_apply:
        assert results_dict[rule]["status"]
    rules_not_apply = [rule for rule in key_list if rule not in rules_apply]
    for rule in rules_not_apply:
        assert not results_dict[rule]["status"]


def test_acmg_RAD51D_intron_variant():
    """
    Test general ACMG guidelines
    """
    path_variant = paths.TEST / "test_variants" / "RAD51D_intron_variant.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    _, results = classify(path_config, var_str)
    results_dict = json.loads(results)
    key_list = [key for key in results_dict.keys()]
    rules_apply = ["PP3_splicing", "PM2"]
    for rule in rules_apply:
        assert results_dict[rule]["status"]
    rules_not_apply = [rule for rule in key_list if rule not in rules_apply]
    for rule in rules_not_apply:
        assert not results_dict[rule]["status"]


def test_acmg_RAD51D_start_loss():
    """
    Test general ACMG guidelines
    """
    path_variant = paths.TEST / "test_variants" / "RAD51D_start_lost.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    _, results = classify(path_config, var_str)
    results_dict = json.loads(results)
    key_list = [key for key in results_dict.keys()]
    rules_apply = ["PVS1_protein", "PM2", "BP4_splicing"]
    for rule in rules_apply:
        assert results_dict[rule]["status"]
    rules_not_apply = [rule for rule in key_list if rule not in rules_apply]
    for rule in rules_not_apply:
        assert not results_dict[rule]["status"]
