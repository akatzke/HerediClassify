#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info
from variant import PopulationDatabases


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
        if flossies.count is None:
            raise ValueError(
                f"The FLOSSIES allele count is None. Please check variant import."
            )
        if flossies.count >= threshold_bs2:
            comment = f"The variant occures {flossies.count} in FLOSSIES (threshold: {threshold_bs2})."
            result = True
        else:
            comment = f"The variant occures {flossies.count} in FLOSSIES (threshold: {threshold_bs2})."
            result = False
        return RuleResult(
            "BS2",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Bs2_with_supporting(abstract_rule):
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
                class_info.THRESHOLD_BS2_SUPPORTING,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        flossies: PopulationDatabases,
        threshold_bs2: int,
        threshold_bs2_supporting: int,
    ) -> RuleResult:
        if flossies.count is None:
            raise ValueError(
                f"The FLOSSIES allele count is None. Please check variant import."
            )
        if flossies.count >= threshold_bs2:
            comment = f"The variant occures {flossies.count} in FLOSSIES (threshold: {threshold_bs2})."
            strength = evidence_strength.STRONG
            result = True
        elif flossies.count >= threshold_bs2_supporting:
            comment = f"The variant occures {flossies.count} in FLOSSIES meeting threshold for supporting evidence strength (threshold: {threshold_bs2_supporting})."
            strength = evidence_strength.SUPPORTING
            result = True
        else:
            comment = f"The variant occures {flossies.count} in FLOSSIES."
            strength = evidence_strength.STRONG
            result = False
        return RuleResult(
            "BS2",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )
