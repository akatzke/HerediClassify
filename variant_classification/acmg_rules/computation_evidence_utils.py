#!/usr/bin/env python3

from typing import Optional
from enum import Enum
from dataclasses import dataclass

from acmg_rules.utils import evidence_strength


class THRESHOLD_DIRECTION(Enum):
    GREATER = "greater"
    GREATER_THAN_OR_EQUAL = "greater_than_or_equal"
    LESS = "less"
    LESS_THAN_OR_EQUAL = "less_than_or_equal"

    @classmethod
    def list(cls):
        """
        Return lisst of all possible Enum values
        """
        return list(map(lambda c: c.value, cls))


@dataclass
class Threshold:
    name: str
    direction: THRESHOLD_DIRECTION
    thresholds: list[float]
    strengths: list[evidence_strength]


def assess_thresholds(
    thresholds: Threshold, prediction_value: Optional[float]
) -> Optional[int]:
    if prediction_value is None:
        return None
    count = 0
    for single_threshold in thresholds.thresholds:
        result = asses_threshold(
            single_threshold, thresholds.direction, prediction_value
        )
        if result is None:
            raise ValueError(
                f"No (valid) direction for the assessment of the prediction tool {thresholds.name} was provided."
            )
        if result:
            count += 1
    return count


# def assess_prediction_tool(
def asses_threshold(
    threshold: float, threshold_direction: THRESHOLD_DIRECTION, prediction_value: float
) -> Optional[bool]:
    """
    Assess prediction result
    """
    if threshold_direction.value == THRESHOLD_DIRECTION.GREATER.value:
        if prediction_value > threshold:
            return True
        else:
            return False
    if threshold_direction.value == THRESHOLD_DIRECTION.GREATER_THAN_OR_EQUAL.value:
        if prediction_value >= threshold:
            return True
        else:
            return False
    elif threshold_direction.value == THRESHOLD_DIRECTION.LESS.value:
        if prediction_value < threshold:
            return True
        else:
            return False
    elif threshold_direction.value == THRESHOLD_DIRECTION.LESS_THAN_OR_EQUAL.value:
        if prediction_value <= threshold:
            return True
        else:
            return False
    else:
        return None


def get_evidence_strength_from_prediction_count(
    strength_order: list[evidence_strength], count: int
) -> evidence_strength:
    """
    Depending on number of thresholds met, return evidence strength
    """
    return strength_order[count - 1]
