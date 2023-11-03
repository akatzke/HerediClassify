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


class pm2(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.VARIANT_GNOMAD,
                info.classification_information.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases, threshold_pm2: float
    ) -> RuleResult:
        if gnomad.frequency > threshold_pm2:
            comment = f"Variant occures with {gnomad.frequency} in {gnomad.name}."
            result = False
        else:
            comment = f"Variant occurs with {gnomad.frequency} in {gnomad.name}."
            result = True
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
