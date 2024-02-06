#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Classification_Info, Info
from variant import PopulationDatabases


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
            (class_info.VARIANT_HOTSPOT_ANNOTATION,),
        )

    @classmethod
    def assess_rule(cls, variant_in_hotspot: bool) -> RuleResult:
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


class Pm1_tp53(abstract_rule):
    """
    PM1: Variant located in mutational hot spot or citical protein region
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_HOTSPOT_ANNOTATION,
                class_info.VARIANT_CANCERHOTSPOTS,
                class_info.THRESHOLD_CANCERHOTSPOTS_AC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant_in_hotspot: bool,
        cancerhotspots: PopulationDatabases,
        threshold_cancerhotspots_ac: float,
    ) -> RuleResult:
        if variant_in_hotspot:
            comment = f"Variant in defined mutational hotspot."
            result = True
        elif cancerhotspots.count is None:
            comment = "Variant not located in defined mutational hotspot and no Cancer Hotspots entry database for the variant."
            result = False
        elif cancerhotspots.count >= threshold_cancerhotspots_ac:
            comment = f"Variant occurs in Cancer Hotspots {cancerhotspots.count} times."
            result = True
        else:
            comment = f"Variant not in mutational hotspot and occurence of variant in Cancer Hotspots ({cancerhotspots.count}) does not meet threshold ({threshold_cancerhotspots_ac})."
            result = False
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
