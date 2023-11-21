#!/usr/bin/env python3

from typing import Optional
from enum import Enum
from dataclasses import dataclass

from variant_classification.acmg_rules.utils import evidence_strength


class THRESHOLD_DIRECTION(Enum):
    HIGHER = "higher"
    LOWER = "lower"


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


def assess_prediction_tool(threshold: Threshold, prediction_value: float) -> bool:
    """
    Assess prediction result
    """
    if threshold.direction.value == THRESHOLD_DIRECTION.HIGHER.value:
        if prediction_value > threshold.threshold:
            return True
        else:
            return False
    elif threshold.direction.value == THRESHOLD_DIRECTION.LOWER.value:
        if prediction_value < threshold.threshold:
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
    if threshold.direction.value == THRESHOLD_DIRECTION.HIGHER.value:
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
    elif threshold.direction.value == THRESHOLD_DIRECTION.LOWER.value:
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
