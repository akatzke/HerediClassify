#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    abstract_rule,
    rule_type,
)
from acmg_rules.computation_evidence_utils import (
    Threshold,
    assess_prediction_tool,
)
from information import Classification_Info, Info


class Pp3_protein(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        try:
            prediction_value = prediction_dict[threshold.name]
            prediction = assess_prediction_tool(threshold, prediction_value)
        except KeyError:
            prediction = None
        if prediction is None:
            comment = f"No score was provided for {threshold.name}"
            result = False
        elif prediction:
            comment = f"Variant is predicted to be pathogenic by {threshold.name}."
            result = True
        else:
            comment = f"Variant is not predicted to be pathogenic by {threshold.name}."
            result = False
        return RuleResult(
            "PP3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Pp3_splicing(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        try:
            prediction_value = prediction_dict[threshold.name]
            prediction = assess_prediction_tool(threshold, prediction_value)
        except KeyError:
            prediction = None
        if prediction is None:
            comment = f"No score was provided for {threshold.name}"
            result = False
        elif prediction:
            comment = (
                f"Variant is predicted to have a splice effect by {threshold.name}."
            )
            result = True
        else:
            comment = (
                f"Variant is not predicted to have a splice effect by {threshold.name}."
            )
            result = False
        return RuleResult(
            "PP3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
