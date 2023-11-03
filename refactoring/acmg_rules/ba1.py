#!/usr/bin/env python3

from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
import refactoring.information as info
from refactoring.variant import PopulationDatabases


class ba1(abstract_rule):
    """
    BA1: High frequency of variant in healthy population (e.g. gnomAD)
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.VARIANT_GNOMAD,
                info.classification_information.THRESHOLD_BA1,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases, threshold_ba1: float
    ) -> RuleResult:
        if gnomad.frequency > threshold_ba1:
            comment = f"Variant occures with {gnomad.frequency} in {gnomad.name}."
            result = True
        else:
            comment = f"Variant occurs with {gnomad.frequency} in {gnomad.name}."
            result = False
        return RuleResult(
            "BA1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STAND_ALONE,
            comment,
        )
