#!/usr/bin/env python3

from typing import Protocol
from enum import Enum

from refactoring.variant_annotate import Variant_annotated
from refactoring.acmg_rules import *


class classification_information(Enum):
    ANNOTATED_TRANSCRIPT_LIST = "annotated_transcript_list"
    VARIANT_CLINVAR = "variant_clinvar"
    VARIANT_HOTSPOT = "variant_hotspot"
    VARIANT_ANNOT = "variant_annot"
    THRESHOLD_PATHOGENICITY_PREDICTION = "threshold_pathogenicity_prediction"
    THRESHOLD_SPLICING_PREDICTION = "threshold_splicing_prediction"
    THRESHOLD_PM2 = "threshold_pm2"
    THRESHOLD_BA1 = "threshold_ba1"
    THRESHOLD_BS1 = "threshold_ps1"
    THRESHOLD_BS2 = "threshold_bs2"
    THRESHOLD_DIFF_LEN_PROT_PERCENT = "threshold_diff_len_prot_percent"


def from_rule_list_get_annotations_needed(rule_list: list[str]) -> list[classification_information]:
    annotations_needed = []
    RULE_DICTIONARY = {"pvs1" : pvs1, "pm1" : pm1, "pm2" : pm2, "pm4" : pm4, "pm5": pm5, "pp3": bp4_pp3, "ba1" : ba1, "bs1": bs1, "bp3" : bp3, "bp4": bp4_pp3, "bp7": bp7}
    for rule in rule_list:
        try:
            rule_class = RULE_DICTIONARY[rule.lower()]
        except:
            print(f"{rule.upper()} not valid rule. Valid rules are {RULE_DICTIONARY.keys()}")
            break
        annotations_needed.append(rule_class.arguments)
    annotations_list = [item for sublist in annotations_needed for item in sublist]
    return annotations_list


def from_annotations_list_get_annotation_function(annotation_list = list[classification_information]) -> list[str]:
    set_annot = set(annotation_list)
    ANNOTATION_DICTIONARY = {classification_information.ANNOTATED_TRANSCRIPT_LIST: }
    return


class Annotation(Protocol):
    """
    Function that annotates a Variant_annotated object
    """

    def __call__(self, variant: Variant_annotated) -> Variant_annotated:
        """
        Return a further annotated Variant_annotated object
        """
        ...
