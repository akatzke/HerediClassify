#!/usr/bin/env python3

import pathlib
import os
import logging

mylogger = logging.getLogger()

os.chdir("/home/katzkean/variant_classification")

from refactoring.config_annotation import perform_annotation
from test.example_variant_splicing_GRCh38 import (
    create_test_variant,
    create_test_variant_no_flossies,
    create_test_variant_no_hotspot,
)


def test_perform_complete_annotation():
    test_var = create_test_variant()
    path_config = pathlib.Path("/home/katzkean/variant_classification/config.yaml")
    rules_list = perform_annotation(path_config, test_var)


def test_perform_annotation_no_path():
    test_var = create_test_variant()
    path_config = pathlib.Path(
        "/home/katzkean/variant_classification/test/config_no_path.yaml"
    )
    rules_list = perform_annotation(path_config, test_var)


def test_perform_annotation_no_prediction():
    test_var = create_test_variant()
    path_config = pathlib.Path(
        "/home/katzkean/variant_classification/test/config_no_prediction.yaml"
    )
    rule_list = perform_annotation(path_config, test_var)


def test_perform_annotation_no_threshold():
    test_var = create_test_variant()
    path_config = pathlib.Path(
        "/home/katzkean/variant_classification/test/config_no_threshold.yaml"
    )
    rule_list = perform_annotation(path_config, test_var)


def test_perform_annotation_wrong_path():
    test_var = create_test_variant()
    path_config = pathlib.Path(
        "/home/katzkean/variant_classification/test/config_wrong_path.yaml"
    )
    rule_list = perform_annotation(path_config, test_var)


def test_perform_annotation_no_flossies():
    test_var = create_test_variant_no_flossies()
    path_config = pathlib.Path("/home/katzkean/variant_classification/config.yaml")
    rules_list = perform_annotation(path_config, test_var)


def test_perform_annotation_no_hotspot():
    test_var = create_test_variant_no_hotspot()
    path_config = pathlib.Path("/home/katzkean/variant_classification/config.yaml")
    rules_list = perform_annotation(path_config, test_var)


# def test_one_threshold_pathogenic() -> None:
#    assert assess_one_threshold(data=1.2, threshold=1.1) == Prediction_result.PATHOGENIC
