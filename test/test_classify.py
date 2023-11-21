#!/usr/bin/env python3

import pathlib

from variant_classification.classify import classify
from test.test_import_variant import create_json_string_from_variant


def test_classify():
    path_variant = pathlib.Path("./API/example_input.json")
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path("./config.yaml")
    results = classify(path_config, variant_str)
