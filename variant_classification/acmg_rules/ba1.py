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
from variant import PopulationDatabases_gnomAD


class Ba1(abstract_rule):
    """
    BA1: High frequency of variant in healthy population (e.g. gnomAD)
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD,
                class_info.THRESHOLD_BA1,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_ba1: float
    ) -> RuleResult:
        if gnomad.popmax_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.popmax_frequency > threshold_ba1:
            comment = f"Variant occures with {gnomad.popmax_frequency} in gnomAD subpopulation {gnomad.popmax}."
            result = True
            if gnomad.popmax == "None":
                comment = f"Variant does not occur in gnomAD, allele frequency in gnomAd is assumed to be 0."
            if gnomad.popmax == "ALL":
                comment = f"Variant has no popmax. Variant occurs with {gnomad.popmax_frequency} in gnomAD."
        else:
            comment = f"Variant occures with {gnomad.popmax_frequency} in gnomAD subpopulation {gnomad.popmax}."
            result = False
        return RuleResult(
            "BA1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STAND_ALONE,
            comment,
        )


class Ba1_with_absolute(abstract_rule):
    """
    BA1: High frequency of variant in healthy population (e.g. gnomAD)
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD,
                class_info.THRESHOLD_BA1,
                class_info.THRESHOLD_BA1_ABSOLUTE,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        gnomad: PopulationDatabases_gnomAD,
        threshold_ba1: float,
        threshold_ba1_absolute: int,
    ) -> RuleResult:
        if gnomad.popmax_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif (
            gnomad.popmax_frequency > threshold_ba1
            and gnomad.popmax_allele_count >= threshold_ba1_absolute
        ):
            comment = f"Variant occures with a frequeny of {gnomad.popmax_frequency} and a total of {gnomad.popmax_allele_count} times in gnomAD subpopulation {gnomad.popmax}."
            result = True
            if gnomad.popmax == "None":
                comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
            if gnomad.popmax == "ALL":
                comment = f"Variant has no popmax. Variant occurs with {gnomad.popmax_frequency} and a total of {gnomad.popmax_allele_count} times in gnomAD."
        else:
            comment = f"Variant occures with {gnomad.popmax_frequency} and a total of {gnomad.popmax_allele_count} times in gnomAD subpopulation {gnomad.popmax}."
            result = False
        return RuleResult(
            "BA1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STAND_ALONE,
            comment,
        )
