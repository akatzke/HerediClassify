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
from var_type import VARTYPE_GROUPS
from variant import TranscriptInfo, VariantInfo


class Pp2(abstract_rule):
    """
    PP2: Missense variant in a gene that has a low rate of benign missense variation and where missense variants are a common mechanism of disease.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT, class_info.TRANSCRIPT),
        )

    @classmethod
    def assess_rule(
        cls, variant: VariantInfo, transcripts: list[TranscriptInfo]
    ) -> RuleResult:
        # In case one disase variant transcripts is defined, use type of variant in that transcript
        # Otherwise use all variant types defined for variant
        if len(transcripts) == 1:
            variant_types = transcripts[0].var_type
        else:
            variant_types = variant.var_type
        if any(var_type in VARTYPE_GROUPS.MISSENSE.value for var_type in variant_types):
            result = True
            comment = f"Missense variants in {variant.gene_name} are a known mechanism of disease."
        else:
            result = False
            comment = f"PP2 does not apply to this variant, as PP2 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
        return RuleResult(
            "PP2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
