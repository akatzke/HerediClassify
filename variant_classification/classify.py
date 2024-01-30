#!/usr/bin/env python3

import pathlib
import argparse

from ensembl import ensembl
from load_config import load_config, get_gene_specific_config
from load_variant import load_variant
from check_disease_relevant_transcript import check_disease_relevant_transcript
from information import Classification_Info
from config_annotation import (
    get_annotations_needed_from_rules,
    get_annotation_functions,
    get_unique_annotations_needed,
    execute_annotation,
    remove_rules_with_missing_annotation,
    apply_rules,
)
from create_output import create_output, create_rules_dict
from check_incompatible_rules import check_incompatible_rules

from os import path
import json
import sys


def classify(config_path: pathlib.Path, variant_str: str):
    """
    Perform classification
    """
    config = load_config(config_path)
    variant = load_variant(variant_str)
    configuration_name, final_config = get_gene_specific_config(
        config, variant.variant_info.gene_name
    )
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    class_info = Classification_Info()
    annotations_needed_by_rules = get_annotations_needed_from_rules(
        final_config["rules"], class_info
    )
    annotations_needed = get_unique_annotations_needed(annotations_needed_by_rules)
    annotations_set_up = get_annotation_functions(
        annotations_needed, variant_disease_relevant, final_config, class_info
    )
    annotation = execute_annotation(annotations_set_up)
    annotations_needed_by_rules_filtered = remove_rules_with_missing_annotation(
        annotations_needed_by_rules
    )
    rule_results = apply_rules(annotations_needed_by_rules_filtered)
    rule_dict = create_rules_dict(rule_results)
    rule_results_checked = check_incompatible_rules(
        rule_dict, configuration_name, final_config["rules"]
    )
    out_result = create_output(rule_results_checked)
    ensembl.clear_cache()
    return configuration_name, out_result


if __name__ == "__main__":
    # define CLI arguments
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
    parser.add_argument(
        "-o",
        "--output",
        default="",
        help="path to output file",
        type=str,
    )
    # read passed CLI arguments
    args = parser.parse_args()

    # check if arguments were given
    if args.input == "":
        input_file = sys.stdin
    if args.config == "":
        raise ValueError("No config file provided.")

    # Execute classification
    path_config = pathlib.Path(args.config)
    input = args.input
    if path.exists(input):
        with open(input) as infile:
            input = infile.read()

    config_name, result = classify(path_config, input)

    # write classification to sout or to file
    if args.output != "":
        sys.stdout = open(args.output, "w")  # overwrite print with sout
    result_json = json.loads(result)
    final_result = json.dumps(result_json, indent=4)
    print(final_result)
