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
from variant import PopulationDatabases_gnomAD, VariantInfo


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
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_pm2: float
    ) -> RuleResult:
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency <= threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no popmax. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
                result = False
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )


class Pm2_supporting(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_pm2: float
    ) -> RuleResult:
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency <= threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no popmax. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Pm2_supporting_faf(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Pm2_supporting.assess_rule,
            (
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_PM2,
            ),
        )


class Pm2_supporting_less(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls, gnomad: PopulationDatabases_gnomAD, threshold_pm2: float
    ) -> RuleResult:
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency < threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no popmax. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no entry in subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Pm2_supporting_less_faf(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Pm2_supporting_less.assess_rule,
            (
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_PM2,
            ),
        )


class Pm2_supporting_no_ins_del_indel(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    PM2 does not apply to insertions, deletions or indels
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_pm2: float,
    ) -> RuleResult:
        if len(variant.var_ref) != 1 or len(variant.var_obs) != 1:
            comment = f"PM2 does not apply to insertions, deletions or delins."
            result = False
        elif gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency <= threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no popmax. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
            if gnomad.subpopulation == "None":
                comment = f"Variant is absent from gnomAD."
            if gnomad.subpopulation == "ALL":
                comment = f"Variant has no entry in subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Pm2_supporting_no_ins_del_indel_faf(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    PM2 does not apply to insertions, deletions or indels
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Pm2_supporting_no_ins_del_indel.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_PM2,
            ),
        )
