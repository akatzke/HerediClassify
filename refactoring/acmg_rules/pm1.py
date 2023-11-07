#!/usr/bin/env python3

from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from refactoring.information import Classification_Info, Info
from refactoring.variant import AffectedRegion


class Pm1(abstract_rule):
    """
    PM1: Variant located in mutational hot spot or citical protein region
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_HOTSPOT,),
        )

    @classmethod
    def assess_rule(cls, variant_in_hotspot: AffectedRegion) -> RuleResult:
        if variant_in_hotspot:
            comment = f"Variant in mutational hotspot."
            result = True
        else:
            comment = "Variant not in mutational hotspot"
            result = False
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
