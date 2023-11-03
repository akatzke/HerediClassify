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


class bs1(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.VARIANT_GNOMAD,
                info.classification_information.THRESHOLD_BS1,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases, threshold_bs1: float
    ) -> RuleResult:
        if gnomad.frequency > threshold_bs1:
            comment = f"Variant occures with {gnomad.frequency} in {gnomad.name}."
            result = True
        else:
            comment = f"Variant occurs with {gnomad.frequency} in {gnomad.name}."
            result = False
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )
