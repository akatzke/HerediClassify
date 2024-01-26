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
from var_type import VARTYPE, VARTYPE_GROUPS
from variant import VariantInfo
from acmg_rules.computation_evidence_utils import Threshold, assess_prediction_tool
from var_type import VARTYPE


class Bp1_annotation_cold_spot_strong(abstract_rule):
    """
    BP1: Missense variant in a gene that has a low rate of benign missense variation and where missense variants are a common mechanism of disease.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_COLDSPOT_ANNOTATION,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        coldspot: bool,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        prediction_benign = assess_prediction_tool(threshold, prediction_value)
        variant_types = [
            VARTYPE.SYNONYMOUS_VARIANT,
            VARTYPE.MISSENSE_VARIANT,
            VARTYPE.INFRAME_DELETION,
            VARTYPE.INFRAME_INSERTION,
        ]
        if (
            coldspot
            and prediction_benign
            and any(var_type in variant_types for var_type in variant.var_type)
        ):
            result = True
            comment = (
                f"Variant is located in coldspot region as defined in annotation file."
            )
        elif not coldspot:
            result = False
            comment = f"Variant is not located in coldspot region as defined in annotation file."
        elif prediction_benign is None:
            result = False
            comment = f"No splicing prediction is available for variant. Therefore, BP1 does not apply."
        elif not prediction_benign:
            result = False
            comment = f"Variant is not predicted to not affect splicing. Therefore, BP1 does not apply."
        elif not any(var_type in variant_types for var_type in variant.var_type):
            result = False
            comment = f"BP1 does not apply to this variant, as BP1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
        else:
            result = False
            comment = f"Variant is not located in coldspot region as defined in annotation file."
        return RuleResult(
            "BP1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Bp1(abstract_rule):
    """
    BP1: Missense variant in a gene for which primarily truncating variants are known to cause disease.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT,),
        )

    @classmethod
    def assess_rule(cls, variant: VariantInfo) -> RuleResult:
        if any(
            var_type in VARTYPE_GROUPS.MISSENSE.value for var_type in variant.var_type
        ):
            result = True
            comment = f"Missense variants in {variant.gene_name} where primarily truncating variants are known to cause disease."
        else:
            result = False
            comment = f"BP1 does not apply to this variant, as BP1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
        return RuleResult(
            "BP1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
