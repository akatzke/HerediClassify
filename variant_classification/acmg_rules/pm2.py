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


class Pm2(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD,
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_pm2: float
    ) -> RuleResult:
        if gnomad.popmax_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.popmax_frequency > threshold_pm2:
            comment = f"Variant occures with {gnomad.popmax_frequency} in GnomAD subpopulation {gnomad.popmax}."
            result = False
        else:
            comment = f"Variant occures with {gnomad.popmax_frequency} in GnomAD subpopulation {gnomad.popmax}."
            result = True
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
