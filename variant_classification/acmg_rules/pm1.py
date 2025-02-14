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
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds


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


class Pm1_supporting(abstract_rule):
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
            evidence_strength.SUPPORTING,
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
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant_in_hotspot: bool,
        cancerhotspots: PopulationDatabases,
        threshold_cancerhotspots_ac: float,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if variant_in_hotspot:
            comment = f"Variant in defined mutational hotspot."
            result = True
        elif cancerhotspots.count is None:
            comment = f"Variant not located in defined mutational hotspot and no Cancer Hotspots entry database for the variant."
            result = False
        elif cancerhotspots.count >= threshold_cancerhotspots_ac:
            comment = f"Variant occurs in Cancer Hotspots {cancerhotspots.count} times."
            result = True
        else:
            comment = f"Variant not in mutational hotspot and occurence of variant in Cancer Hotspots ({cancerhotspots.count}) does not meet threshold ({threshold_cancerhotspots_ac})."
            result = False
        if num_thresholds_met is None:
            comment = (
                "ATTENTION: No splicing prediction available for this variant. "
                + comment
            )
        elif num_thresholds_met == 0:
            comment = f"Variant is not predicted to not affect splicing, therefore PM1 does not apply."
            result = False
        else:
            comment = (
                comment
                + f" {threshold.name} predicts no effect on splicing for the variant (threshold: {threshold.thresholds[num_thresholds_met -1]}, value: {prediction_value})."
            )
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
