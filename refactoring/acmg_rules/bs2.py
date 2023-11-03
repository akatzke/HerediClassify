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


class bs2(abstract_rule):
    """
    BS2: Mutation found in healthy individuals
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.VARIANT_FLOSSIES,
                info.classification_information.THRESHOLD_BS1,
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
