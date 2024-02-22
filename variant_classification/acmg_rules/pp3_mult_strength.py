#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    abstract_rule,
    rule_type,
)
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from information import Classification_Info, Info


class Pp3_protein_mult_strength(abstract_rule):
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
            num_thresholds_met = assess_thresholds(threshold, prediction_value)
        except KeyError:
            num_thresholds_met = None
        if num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}"
            strength = evidence_strength.SUPPORTING
            result = False
        elif num_thresholds_met == 0:
            comment = f"Variant is not predicted to be pathogenic by {threshold.name}."
            strength = evidence_strength.SUPPORTING
            result = False
        else:
            result = True
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to be pathogenic by {threshold.name} with evidence strength {strength.value}."
        return RuleResult(
            "PP3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )


class Pp3_splicing_mult_strength(abstract_rule):
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
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}"
            result = False
            strength = evidence_strength.SUPPORTING
        elif num_thresholds_met == 0:
            comment = (
                f"Variant is not predicted to have a splice effect by {threshold.name}."
            )
            result = False
            strength = evidence_strength.SUPPORTING
        else:
            result = True
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to have a splice effect by {threshold.name} with evidence strength {strength.value}."
        return RuleResult(
            "PP3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
