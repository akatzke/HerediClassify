#!/usr/bin/env python3
import yaml
from typing import Literal
from dacite import

from refactoring.rules import STRENGTH_TYPE, assess_pp3

STRENGTH_TYPE = Literal["supporting", "moderate", "strong", "very strong"]

path_config = "/home/katzkean/variant_classification/config.yaml"

def check_dictionary_for_function(rules)

@dataclass
class Configuration():
    rules: dict[str, dict]
    prediction_tools_threshold: dict[str, dict]
    allele_frequency_threshold: dict[str, int]


    def load_config(path_yaml: str) -> Configuration:
        with open(path_yaml) as config_file:
            config = yaml.safe_load(config_file)
         return config

@dataclass
class Prediction_too_threshold():
    predictions: dict[str, int]

@dataclass
class Allele_frequency_thresholds():
    thresholds: dict[str, int]
