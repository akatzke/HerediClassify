#!/usr/bin/env python3

from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
import refactoring.information as info
from refactoring.clinvar_utils import CLINVAR_TYPE, ClinVar


class ps1_protein(abstract_rule):
    """
    PS1: Position is classified as pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (info.classification_information.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_result: dict[CLINVAR_TYPE, ClinVar]) -> RuleResult:
        clinvar_same_aa = clinvar_result[CLINVAR_TYPE.SAME_AA_CHANGE]
        if clinvar_same_aa.pathogenic:
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}."
            result = True
        else:
            comment = "No matches found for variant."
            result = False
        return RuleResult(
            "PS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )


class ps1_splicing(abstract_rule):
    """
    PS1 for splicing: Splice variant in same position has been show to be pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (info.classification_information.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_result: dict[CLINVAR_TYPE, ClinVar]) -> RuleResult:
        clinvar_same_nucleotide = clinvar_result[CLINVAR_TYPE.SAME_NUCLEOTIDE]
        if clinvar_same_nucleotide.pathogenic:
            comment = f"The following ClinVar entries show splice variants at the same nucleotide position to be pathogenic: {clinvar_same_nucleotide.ids}."
            result = True
        else:
            comment = "No matches found for variant."
            result = False
        return RuleResult(
            "PS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )
