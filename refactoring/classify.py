#!/usr/bin/env python3

import pathlib
import argparse

from refactoring.variant import Variant
from refactoring.config_annotation import perform_annotation

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


def classify(config_path: pathlib.Path, variant: pathlib.Path) -> :
    """
    Perform classification
    """
    pass


if __name__ == "__main__":
    path_config = pathlib.Path(args.config)
    path_input = pathlib.Path(args.input)
    classify(path_config, path_input)
