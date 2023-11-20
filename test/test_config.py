#!/usr/bin/env python3

import pathlib
import logging

mylogger = logging.getLogger()

from refactoring.classify import classify
from test.test_import_variant import create_json_string_from_variant


def test_classify_no_path():
    path_variant = pathlib.Path("./API/example_input.json")
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path("./test/config_no_path.yaml")
    rules_list = classify(path_config, variant_str)


def test_classify_no_prediction():
    path_variant = pathlib.Path("./API/example_input.json")
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path("./test/config_no_prediction.yaml")
    rule_list = classify(path_config, variant_str)


def test_classify_no_threshold():
    path_variant = pathlib.Path("./API/example_input.json")
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path("./test/config_no_threshold.yaml")
    rule_list = classify(path_config, variant_str)


def test_classify_wrong_path():
    path_variant = pathlib.Path("./API/example_input.json")
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path("./test/config_wrong_path.yaml")
    rule_list = classify(path_config, variant_str)


def test_classify_no_flossies():
    path_variant = pathlib.Path("./API/example_input.json")
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path("./config.yaml")
    rules_list = classify(path_config, variant_str)


def test_classify_no_hotspot():
    path_variant = pathlib.Path("./API/example_input.json")
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path("./config.yaml")
    rules_list = classify(path_config, variant_str)


# def test_one_threshold_pathogenic() -> None:
#    assert assess_one_threshold(data=1.2, threshold=1.1) == Prediction_result.PATHOGENIC
