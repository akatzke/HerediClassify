#!/usr/bin/env python3

import pathlib
import json

from refactoring.load_variant import load_variant


def create_json_string_from_variant(path: pathlib.Path) -> str:
    """
    Load json file and return json string
    """
    with open(path) as f:
        variant = json.load(f)
    json_string = json.dumps(variant)
    return json_string


def test_variant_import():
    path_variant = pathlib.Path(
        "/home/katzkean/variant_classification/API/example_input.json"
    )
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)


def test_variant_import_missing_flossies():
    path_variant = pathlib.Path(
        "/home/katzkean/variant_classification/test/test_variants/test_var_no_flossies.json"
    )
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)


def test_variant_import_missing_cancer_hotspot():
    path_variant = pathlib.Path(
        "/home/katzkean/variant_classification/test/test_variants/test_var_no_cancer_hotspot.json"
    )
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)


def test_variant_import_missing_chr():
    path_variant = pathlib.Path(
        "/home/katzkean/variant_classification/test/test_variants/test_var_no_chr.json"
    )
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
