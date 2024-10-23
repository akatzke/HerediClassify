#!/usr/bin/env python3

from typing import Callable
import pandas as pd
from functools import lru_cache


def get_final_classification_from_possible_classes(possible_class: list[int]) -> int:
    """
    Get final classifcation from list of possible classes
    """
    classes_set = set(possible_class)
    if not classes_set:
        return 3
    elif len(classes_set) == 1:
        return classes_set.pop()
    elif (1 in classes_set or 2 in classes_set) and (
        4 in classes_set or 5 in classes_set
    ):
        return 3
    elif 1 in classes_set and 2 in classes_set:
        return 1
    elif 5 in classes_set and 4 in classes_set:
        return 5


def get_classifications_from_rule_combinations(
    rule_combinations_dict: dict[int, list[Callable]], rule_results: pd.DataFrame
) -> list[int]:
    """
    For every rule recommendation in dictionary check if rule combination applies and add class to output list
    """
    possible_classes = []
    for classification, class_rule_combinations in rule_combinations_dict.items():
        for check_rule_combination in class_rule_combinations:
            if check_rule_combination(rule_results):
                possible_classes.append(classification)
    return possible_classes


@lru_cache
def get_counts(rules: dict[str, int]) -> dict[str, int]:
    df = pd.DataFrame(rules)
    counts = df.groupby("evidence-strength").count()
    return counts.to_dict()


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
    def fun(rules: dict[str, int]) -> bool:
        counts = get_counts(rules)
        prelim_results = []
        if min_benign_stand_alone is not None:
            prelim_results.append(
                counts["benign_stand_alone"] >= min_benign_stand_alone
            )
        if min_benign_strong is not None:
            prelim_results.append(counts["benign_strong"] >= min_benign_strong)
        if min_benign_moderate is not None:
            prelim_results.append(counts["benign_moderate"] >= min_benign_moderate)
        if min_benign_supporting is not None:
            prelim_results.append(counts["benign_supporting"] >= min_benign_supporting)
        if min_pathogenic_very_strong is not None:
            prelim_results.append(
                counts["pathogenic_very_strong"] >= min_pathogenic_very_strong
            )
        if min_pathogenic_strong is not None:
            prelim_results.append(counts["pathogenic_strong"] >= min_pathogenic_strong)
        if min_pathogenic_moderate is not None:
            prelim_results.append(
                counts["pathogenic_moderate"] >= min_pathogenic_moderate
            )
        if min_pathogenic_supporting is not None:
            prelim_results.append(
                counts["pathogenic_supporting"] >= min_pathogenic_supporting
            )
        return all(prelim_results)

    return fun
