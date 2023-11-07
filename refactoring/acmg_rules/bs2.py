#!/usr/bin/env python3

from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from refactoring.information import Classification_Info, Info
from refactoring.variant import PopulationDatabases


class Bs2(abstract_rule):
    """
    BS2: Mutation found in healthy individuals
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_FLOSSIES,
                class_info.THRESHOLD_BS2,
            ),
        )

    @classmethod
    def assess_rule(
        cls, flossies: PopulationDatabases, threshold_bs2: float
    ) -> RuleResult:
        if flossies.frequency > threshold_bs2:
            comment = "Something"
            result = True
        else:
            comment = "Something"
            result = False
        return RuleResult(
            "BS2",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )
