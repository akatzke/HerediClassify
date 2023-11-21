#!/usr/bin/env python3

from typing import Callable

from variant_classification.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from variant_classification.information import Classification_Info, Info
from variant_classification.acmg_rules.computation_evidence_utils import (
    assess_prediction_tool,
    Threshold,
)


class Bp4_protein(abstract_rule):
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
        try:
            prediction_value = prediction_dict[threshold.name]
        except KeyError:
            raise KeyError(
                f"For {threshold.name} no prediction value was found in {prediction_dict}"
            )
        prediction = assess_prediction_tool(threshold, prediction_value)
        if prediction:
            comment = "Varinat is predicted to be pathogenic."
            result = True
        else:
            comment = "Varinat is not predicted to be pathogenic."
            result = False
        return RuleResult(
            "BP4",
            rule_type.PROTEIN,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Bp4_splicing(abstract_rule):
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
        try:
            prediction_value = prediction_dict[threshold.name]
        except KeyError:
            raise KeyError(
                f"For {threshold.name} no prediction value was found in {prediction_dict}"
            )
        prediction = assess_prediction_tool(threshold, prediction_value)
        if prediction:
            comment = "Varinat is predicted to be benign."
            result = True
        else:
            comment = "Varinat is not predicted to be benign."
            result = False
        return RuleResult(
            "BP4",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
