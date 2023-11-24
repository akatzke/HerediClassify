#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Classification_Info, Info
from clinvar_utils import ClinVar_Type, ClinVar


class Ps1_protein(abstract_rule):
    """
    PS1: Position is classified as pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_result: dict[ClinVar_Type, ClinVar]) -> RuleResult:
        clinvar_same_aa = clinvar_result[ClinVar_Type.SAME_AA_CHANGE]
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


class Ps1_splicing(abstract_rule):
    """
    PS1 for splicing: Splice variant in same position has been show to be pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_result: dict[ClinVar_Type, ClinVar]) -> RuleResult:
        clinvar_same_nucleotide = clinvar_result[ClinVar_Type.SAME_NUCLEOTIDE]
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
