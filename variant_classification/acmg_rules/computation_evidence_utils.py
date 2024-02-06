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
    threshold: float
    direction: THRESHOLD_DIRECTION


@dataclass
class Threshold_evidence_strength:
    name: str
    direction: THRESHOLD_DIRECTION
    threshold_very_strong: Optional[float] = None
    threshold_strong: Optional[float] = None
    threshold_moderate: Optional[float] = None
    threshold_supporting: Optional[float] = None


def assess_prediction_tool(
    threshold: Threshold, prediction_value: Optional[float]
) -> Optional[bool]:
    """
    Assess prediction result
    """
    if prediction_value is None:
        return None
    if threshold.direction.value == THRESHOLD_DIRECTION.GREATER.value:
        if prediction_value > threshold.threshold:
            return True
        else:
            return False
    if threshold.direction.value == THRESHOLD_DIRECTION.GREATER_THAN_OR_EQUAL.value:
        if prediction_value >= threshold.threshold:
            return True
        else:
            return False
    elif threshold.direction.value == THRESHOLD_DIRECTION.LESS.value:
        if prediction_value < threshold.threshold:
            return True
        else:
            return False
    elif threshold.direction.value == THRESHOLD_DIRECTION.LESS_THAN_OR_EQUAL.value:
        if prediction_value <= threshold.threshold:
            return True
        else:
            return False
    else:
        raise ValueError(
            f"No direction for the assessment of the prediction tool {threshold.name} was provided."
        )


def assess_prediction_tool_diff_evidence_strength(
    threshold: Threshold_evidence_strength, prediction_value: float
) -> tuple[bool, evidence_strength]:
    """
    Return
    """
    if threshold.direction.value == THRESHOLD_DIRECTION.GREATER.value:
        if (
            threshold.threshold_very_strong is not None
            and prediction_value > threshold.threshold_very_strong
        ):
            return (True, evidence_strength.VERY_STRONG)
        elif (
            threshold.threshold_strong is not None
            and prediction_value > threshold.threshold_strong
        ):
            return (True, evidence_strength.STRONG)
        elif (
            threshold.threshold_moderate is not None
            and prediction_value > threshold.threshold_moderate
        ):
            return (True, evidence_strength.MODERATE)
        elif (
            threshold.threshold_supporting is not None
            and prediction_value > threshold.threshold_supporting
        ):
            return (True, evidence_strength.SUPPORTING)
        else:
            return (False, evidence_strength.SUPPORTING)
    elif threshold.direction.value == THRESHOLD_DIRECTION.LESS.value:
        if (
            threshold.threshold_very_strong is not None
            and prediction_value < threshold.threshold_very_strong
        ):
            return (True, evidence_strength.VERY_STRONG)
        elif (
            threshold.threshold_strong is not None
            and prediction_value < threshold.threshold_strong
        ):
            return (True, evidence_strength.STRONG)
        elif (
            threshold.threshold_moderate is not None
            and prediction_value < threshold.threshold_moderate
        ):
            return (True, evidence_strength.MODERATE)
        elif (
            threshold.threshold_supporting is not None
            and prediction_value < threshold.threshold_supporting
        ):
            return (True, evidence_strength.SUPPORTING)
        else:
            return (False, evidence_strength.SUPPORTING)
    else:
        raise ValueError(
            f"No direction for the assessment of the prediction tool {threshold.name} was provided."
        )
