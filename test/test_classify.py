from variant_classification.classify import classify
from test.test_import_variant import create_json_string_from_variant
import test.paths as paths


def test_classify():
    path_variant = paths.API / "example_input.json"
    variant_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    results = classify(path_config, variant_str)


def test_classify_brca1():
    path_variant = paths.TEST / "test_variants" / "test_var_BRCA1_exon8.json"
    variant_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    results = classify(path_config, variant_str)


def test_classify_brca1_2():
    path_variant = paths.TEST / "test_variants" / "test_var_stop_gained_BRCA1_2.json"
    variant_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "config.yaml"
    results = classify(path_config, variant_str)
