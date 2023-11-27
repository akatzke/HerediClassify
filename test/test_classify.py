import pathlib

from os import path
from variant_classification.classify import classify
from test.test_import_variant import create_json_string_from_variant
import test.paths as paths


def test_classify():
    path_variant = pathlib.Path(path.join(paths.API, "example_input.json"))
    variant_str = create_json_string_from_variant(path_variant)
    path_config = pathlib.Path(path.join(paths.ROOT, "config.yaml"))
    results = classify(path_config, variant_str)
