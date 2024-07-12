#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info
from acmg_rules.computation_evidence_utils import (
    assess_thresholds,
    Threshold,
)


class Bp4_protein_mult_strength(abstract_rule):
    """
    BP4: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_PATHOGENICITY_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}."
            result = False
            strength = evidence_strength.SUPPORTING
        elif num_thresholds_met == 0:
            comment = f"Variant is not predicted to be benign by {threshold.name}."
            result = False
            strength = evidence_strength.SUPPORTING
        else:
            result = True
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to be benign by {threshold.name} with evidence strength {strength.value} meeting a threshold of {threshold.thresholds[num_thresholds_met -1]} (value: {prediction_value})."
        return RuleResult(
            "BP4",
            rule_type.PROTEIN,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )


class Bp4_splicing_mult_strength(abstract_rule):
    """
    BP4: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}."
            result = False
            strength = evidence_strength.SUPPORTING
        elif num_thresholds_met == 0:
            comment = f"Variant is not predicted to have no splicing effect by {threshold.name}."
            result = False
            strength = evidence_strength.SUPPORTING
        else:
            result = False
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to have no splicing effect by {threshold.name} with evidence strength {strength.value} meeting a threshold of {threshold.thresholds[num_thresholds_met -1]} (value: {prediction_value})."
        return RuleResult(
            "BP4",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
