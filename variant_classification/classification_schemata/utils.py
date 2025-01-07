#!/usr/bin/env python3

from typing import Callable
import pandas as pd


def get_final_classification_from_possible_classes(possible_class: list[int]) -> int:
    """
    Get final classifcation from list of possible classes
    """
    classes_set = set(possible_class)
    if not classes_set:
        return 3
    if len(classes_set) == 1:
        return classes_set.pop()
    if (1 in classes_set or 2 in classes_set) and (
        4 in classes_set or 5 in classes_set
    ):
        return 3
    if 1 in classes_set and 2 in classes_set:
        return 1
    if 5 in classes_set and 4 in classes_set:
        return 5
    return 3


def get_classifications_from_rule_combinations(
    rule_combinations_dict: dict[int, list[Callable]],
    counts_evidence_strength: dict[str, int],
) -> list[int]:
    """
    For every rule recommendation in dictionary check if rule combination applies and add class to output list
    """
    possible_classes = []
    for classification, class_rule_combinations in rule_combinations_dict.items():
        for check_rule_combination in class_rule_combinations:
            if check_rule_combination(counts_evidence_strength):
                possible_classes.append(classification)
    return possible_classes


def generate_check_specific_rule(rule_name: str) -> Callable[[list], bool]:
    """
    Create function that checks if specific rule applies
    """

    def fun(dict) -> bool:
        applicable_rules = dict.get("applicable_rules", [])
        return rule_name in applicable_rules

    return fun


def generate_count_rule(
    min_benign_stand_alone: int | None = None,
    min_benign_strong: int | None = None,
    min_benign_moderate: int | None = None,
    min_benign_supporting: int | None = None,
    min_pathogenic_very_strong: int | None = None,
    min_pathogenic_strong: int | None = None,
    min_pathogenic_moderate: int | None = None,
    min_pathogenic_supporting: int | None = None,
) -> Callable[[dict[str, int]], int | None]:
    """
    Create function that counts how often rules of a given evidence strength apply
    """

    def fun(counts: dict[str, int]) -> bool:
        prelim_results = []
        if min_benign_stand_alone is not None:
            prelim_results.append(
                counts.get("benign_stand_alone", 0) >= min_benign_stand_alone
            )
        if min_benign_strong is not None:
            prelim_results.append(counts.get("benign_strong", 0) >= min_benign_strong)
        if min_benign_moderate is not None:
            prelim_results.append(
                counts.get("benign_moderate", 0) >= min_benign_moderate
            )
        if min_benign_supporting is not None:
            prelim_results.append(
                counts.get("benign_supporting", 0) >= min_benign_supporting
            )
        if min_pathogenic_very_strong is not None:
            prelim_results.append(
                counts.get("pathogenic_very_strong", 0) >= min_pathogenic_very_strong
            )
        if min_pathogenic_strong is not None:
            prelim_results.append(
                counts.get("pathogenic_strong", 0) >= min_pathogenic_strong
            )
        if min_pathogenic_moderate is not None:
            prelim_results.append(
                counts.get("pathogenic_moderate", 0) >= min_pathogenic_moderate
            )
        if min_pathogenic_supporting is not None:
            prelim_results.append(
                counts.get("pathogenic_supporting", 0) >= min_pathogenic_supporting
            )
        return all(prelim_results)

    return fun
