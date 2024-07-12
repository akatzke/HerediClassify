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
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BA1,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_ba1: float
    ) -> RuleResult:
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency > threshold_ba1:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry in subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "BA1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STAND_ALONE,
            comment,
        )


class Ba1_faf(abstract_rule):
    """
    BA1: High frequency of variant in healthy population (e.g. gnomAD)
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Ba1.assess_rule,
            (
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_BA1,
            ),
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
                class_info.VARIANT_GNOMAD_POPMAX,
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
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif (
            gnomad.subpopulation_frequency > threshold_ba1
            and gnomad.subpopulation_allele_count >= threshold_ba1_absolute
        ):
            comment = f"Variant occures with a frequeny of {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry in subpopulation. Variant occurs with {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD."
        return RuleResult(
            "BA1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STAND_ALONE,
            comment,
        )
