#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Info, Classification_Info
from variant import PopulationDatabases_gnomAD


class Bs1(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BS1,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_bs1: float
    ) -> RuleResult:
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency > threshold_bs1:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Bs1_with_absolute(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BS1,
                class_info.THRESHOLD_BS1_ABSOLUTE,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        gnomad: PopulationDatabases_gnomAD,
        threshold_bs1: float,
        threshold_bs1_absolute: int,
    ) -> RuleResult:
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif (
            gnomad.subpopulation_frequency > threshold_bs1
            and gnomad.subpopulation_allele_count >= threshold_bs1_absolute
        ):
            comment = f"Variant occures with a frequeny of {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            if gnomad.subpopulation == "None":
                comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no entry for subopulation. Variant occurs with {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD."
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Bs1_with_supporting(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BS1,
                class_info.THRESHOLD_BS1_SUPPORTING,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        gnomad: PopulationDatabases_gnomAD,
        threshold_bs1: float,
        threshold_bs1_supporting: float,
    ) -> RuleResult:
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency > threshold_bs1:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            strength = evidence_strength.STRONG
        elif gnomad.subpopulation_frequency > threshold_bs1_supporting:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            strength = evidence_strength.SUPPORTING
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
            strength = evidence_strength.STRONG
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )


class Bs1_with_supporting_faf(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Bs1_with_supporting.assess_rule,
            (
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_BS1,
                class_info.THRESHOLD_BS1_SUPPORTING,
            ),
        )
