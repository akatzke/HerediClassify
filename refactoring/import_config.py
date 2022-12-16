#!/usr/bin/env python3
import yaml
from typing import Literal, Union, Optional
from dataclasses import dataclass
from dacite import from_dict

STRENGTH_TYPE = Literal[
    "supporting", "moderate", "strong", "very strong", "stand_alone"
]

RULE_TYPE = Literal["pathogenic", "benige", "VUS"]

PREDICTION_TYPE = Literal["splicing", "pathogenicity"]

GENE = Literal["BRCA1", "BRCA2", "CHEK2", "BARD1"]


@dataclass
class Prediction_tool_threshold:
    revel_benign: float
    revel_pathogenic: float
    revel: float
    CADD: float
    pyhlop: float
    SpliceAI: float
    MaxEntScan: float
    Hbond: float


@dataclass
class threshold:
    threshold: float
    prediction_tool: str
    prediction_type: PREDICTION_TYPE
    gene: Optional[GENE]


@dataclass
class Allele_frequency_threshold:
    BRCA1: float = 0.1
    BRCA2: float = 0.1
    BARD1: float = 0.1
    BRIP1: float = 0.1
    CDH1: float = 0.1
    CHEK2: float = 0.1
    ATM: float = 0.1
    PALB2: float = 0.1
    PTEN: float = 0.1
    RAD51C: float = 0.1
    RAD51D: float = 0.1
    P53: float = 0.1


@dataclass
class Configuration:
    prediction_tool_threshold: Prediction_tool_threshold
    rules: list[str]
    allele_frequency_threshold: Allele_frequency_threshold


def load_config(path_config: str):
    with open(path_config) as config_file:
        config = yaml.safe_load(config_file)
    return from_dict(Configuration, config)
