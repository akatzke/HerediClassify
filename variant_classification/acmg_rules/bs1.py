#!/usr/bin/env python3

from typing import Callable

from variant_classification.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from variant_classification.information import Info, Classification_Info
from variant_classification.variant import PopulationDatabases_gnomAD


class Bs1(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD,
                class_info.THRESHOLD_BS1,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_bs1: float
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
