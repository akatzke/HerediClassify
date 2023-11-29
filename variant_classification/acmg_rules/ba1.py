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
            comment = f"Variant occures with {gnomad.popmax_frequency} in GnomAD subpopulation {gnomad.popmax}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.popmax_frequency} in GnomAD subpopulation {gnomad.popmax}."
            result = False
        return RuleResult(
            "BA1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STAND_ALONE,
            comment,
        )
