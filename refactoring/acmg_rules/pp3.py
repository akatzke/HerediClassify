#!/usr/bin/env python3

from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    abstract_rule,
    rule_type,
)
from refactoring.acmg_rules.computation_evidence_utils import (
    Threshold,
    assess_prediction_tool,
)
import refactoring.information as info


class pp3_protein(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.VARIANT_PREDICTION,
                info.classification_information.THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC,
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
        except KeyError:
            raise KeyError(
                f"For {threshold.name} no prediction value was found in {prediction_dict}"
            )
        prediction = assess_prediction_tool(threshold, prediction_value)
        if prediction:
            comment = "Variant is predicted to be pathogenic."
            result = True
        else:
            comment = "Variant is not predicted to be pathogenic."
            result = False
        return RuleResult(
            "PP3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class pp3_splicing(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.VARIANT_PREDICTION,
                info.classification_information.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
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
        except KeyError:
            raise KeyError(
                f"For {threshold.name} no prediction value was found in {prediction_dict}"
            )
        prediction = assess_prediction_tool(threshold, prediction_value)
        if prediction:
            comment = "Variant is predicted to be pathogenic."
            result = True
        else:
            comment = "Variant is not predicted to be pathogenic."
            result = False
        return RuleResult(
            "PP3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
