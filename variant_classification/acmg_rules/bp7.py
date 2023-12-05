#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info
from acmg_rules.computation_evidence_utils import (
    assess_prediction_tool,
    Threshold,
)
from variant_classification.variant import VariantInfo
from variant_classification.var_type import VARTYPE


class Bp7(abstract_rule):
    """
    BP7: Silent missense variant is predicted to have effect on splicing
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        if not any(
            var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant.var_type
        ):
            return RuleResult(
                "BP7",
                rule_type.SPLICING,
                evidence_type.BENIGN,
                False,
                evidence_strength.SUPPORTING,
                f"BP7 does not apply to this variant, as PVS1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}.",
            )
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
                f"Variant is predicted to have no splicing effect by {threshold.name}."
            )
            result = True
        else:
            comment = f"Variant is not predicted to have no splicing effect by {threshold.name}."
            result = False
        return RuleResult(
            "BP7",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
