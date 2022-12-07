#!/usr/bin/env python3
import yaml
from typing import Literal
from dataclasses import dataclass
from dacite import from_dict

STRENGTH_TYPE = Literal["supporting", "moderate", "strong", "very strong"]

path_config = "/home/katzkean/variant_classification/config.yaml"

# def check_dictionary_for_function(rules):


@dataclass
class Configuration:
    prediction_tool_threshold: dict[str, dict]
    rules: dict[str, dict]
    # allele_frequency_threshold: dict[str, int]

    def load_config(path_yaml: str):
        with open(path_yaml) as config_file:
            config = yaml.safe_load(config_file)
        test = dacite.from_dict(config)
        return test


@dataclass
class Prediction_too_threshold:
    predictions: dict[str, int]


@dataclass
class Allele_frequency_thresholds:
    thresholds: dict[str, int]


@dataclass
class rules_config:
    rules: dict[str, str]


@dataclass
class A:
    x: str


def test_dacite():
    @dataclass
    class X:
        i: str

    test = {"i": "b"}

    result = from_dict(X, test)
    return result
