#!/usr/bin/env python3

import pathlib
import argparse

from variant_classification.load_config import load_config, get_gene_specific_config
from variant_classification.load_variant import load_variant
from variant_classification.information import Classification_Info
from variant_classification.config_annotation import (
    get_annotations_needed_from_rules,
    get_annotation_functions,
    get_unique_annotations_needed,
    execute_annotation,
    remove_rules_with_missing_annotation,
    apply_rules,
)
from variant_classification.create_output import create_output

parser = argparse = argparse.ArgumentParser()

parser.add_argument(
    "-i", "--input", default="", help="Json string of variant", type=str
)
parser.add_argument(
    "-c",
    "--config",
    default="",
    help="path to configuration",
    type=str,
)
args = parser.parse_args()


def classify(config_path: pathlib.Path, variant_str: str):
    """
    Perform classification
    """
    config = load_config(config_path)
    variant = load_variant(variant_str)
    final_config = get_gene_specific_config(config, variant.variant_info.gene_name)
    class_info = Classification_Info()
    annotations_needed_by_rules = get_annotations_needed_from_rules(
        final_config["rules"], class_info
    )
    annotations_needed = get_unique_annotations_needed(annotations_needed_by_rules)
    annotations_set_up = get_annotation_functions(
        annotations_needed, variant, final_config, class_info
    )
    annotation = execute_annotation(annotations_set_up)
    annotations_needed_by_rules_filtered = remove_rules_with_missing_annotation(
        annotations_needed_by_rules
    )
    rule_results = apply_rules(annotations_needed_by_rules_filtered)
    out_result = create_output(rule_results)
    return out_result


if __name__ == "__main__":
    if args.input == "":
        raise ValueError("No variant json string provided.")
    if args.config == "":
        raise ValueError("No config file provided.")

    path_config = pathlib.Path(args.config)
    classify(path_config, args.input)
