#!/usr/bin/env python3

import pathlib
import argparse
import sys

from refactoring.load_config import load_config, get_gene_specific_config
from refactoring.load_variant import load_variant
from refactoring.config_annotation import (
    get_annotations_from_rule,
    perform_annotation,
    apply_rules,
)

parser = argparse = argparse.ArgumentParser()

parser.add_argument("-i", "--input", default="", help="path to input file", type=str)
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
    annotations_needed = get_annotations_from_rules(final_config, variant)
    annotation = perform_annotation(config, variant)
    rule_result = apply_rules(annotation)
    return rule_result


if __name__ == "__main__":
    if args.output != "":
        sys.stdout = open(args.output, "w")

    if args.input != "":
        input_file = open(args.input, "r")
    else:
        input_file = sys.stdin

    path_config = pathlib.Path(args.config)
    path_input = pathlib.Path(args.input)
    classify(path_config, path_input)
