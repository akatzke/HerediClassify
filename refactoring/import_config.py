#!/usr/bin/env python3
import yaml
from typing import Literal, Union, Optional
from dataclasses import dataclass
from dacite import from_dict

STRENGTH_TYPE = Literal[
    "supporting", "moderate", "strong", "very strong", "stand_alone"
]

RULE_TYPE = Literal["pathogenic", "benige", "VUS"]


@dataclass
class Prediction_tool_threshold:
    revel_benign: Optional[float]
    revel_pathogenic: Optional[float]
    revel: Optional[float]
    CADD: Optional[float]
    pyhlop: Optional[float]
    SpliceAI: Optional[float]
    MaxEntScan: Optional[float]
    Hbond: Optional[float]


@dataclass
class Allele_frequency_threshold:
    BRCA1: Optional[float]
    BRCA2: Optional[float]
    BARD1: Optional[float]
    BRIP1: Optional[float]
    CDH1: Optional[float]
    CHEK2: Optional[float]
    ATM: Optional[float]
    PALB2: Optional[float]
    PTEN: Optional[float]
    RAD51C: Optional[float]
    RAD51D: Optional[float]
    P53: Optional[float]


@dataclass
class Configuration:
    prediction_tool_threshold: Prediction_tool_threshold
    rules: dict[str, list[Union[str, int]]]
    allele_frequency_threshold: Allele_frequency_threshold


def load_config(path_config: str):
    with open(path_config) as config_file:
        config = yaml.safe_load(config_file)
    return from_dict(Configuration, config)
