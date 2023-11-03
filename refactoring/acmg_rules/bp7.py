#!/usr/bin/env python3

from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
import refactoring.information as info
from refactoring.acmg_rules.computation_evidence_utils import (
    assess_prediction_tool,
    Threshold,
)


class bp7(abstract_rule):
    """
    BP7: Silent missense variant is predicted to have effect on splicing
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.VARIANT_PREDICTION,
                info.classification_information.THRESHOLD_SPLICING_PREDICTION_BENIGN,
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
            "BP7",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
