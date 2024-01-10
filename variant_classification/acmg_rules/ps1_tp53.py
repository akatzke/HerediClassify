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
from clinvar_utils import ClinVar_Status, ClinVar_Type, ClinVar
from acmg_rules.computation_evidence_utils import Threshold, assess_prediction_tool
from variant import FunctionalData


class Ps1_protein_tp53(abstract_rule):
    """
    PS1: Position is classified as pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_SPLICEAI_PROTEIN,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
                class_info.SPLICING_ASSAY,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_result: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splice_assay: FunctionalData,
    ) -> RuleResult:
        try:
            prediction_value = prediction_dict[threshold.name]
            prediction = assess_prediction_tool(threshold, prediction_value)
        except KeyError:
            prediction = None
        clinvar_same_aa = clinvar_result[ClinVar_Type.SAME_AA_CHANGE]
        if (
            clinvar_same_aa.pathogenic
            and splice_assay.benign
            and clinvar_same_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}. An RNA shows that the variant does not affect splicing."
            strength = evidence_strength.STRONG
            result = True
        elif (
            clinvar_same_aa.pathogenic
            and not prediction
            and clinvar_same_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}. SpliceAI does not predict an effect on splicing for this variant."
            strength = evidence_strength.MODERATE
            result = True
        elif prediction:
            comment = f"Variant is predicted to affect splicing"
            strength = evidence_strength.STRONG
            result = False
        else:
            comment = "No ClinVar entries found that show the same amino acid change as pathogneic."
            strength = evidence_strength.STRONG
            result = False
        return RuleResult(
            "PS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )


class Ps1_splicing_tp53(abstract_rule):
    """
    PS1 for splicing: Splice variant in same position has been show to be pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_SPLICEAI_SPLICE,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_result: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        try:
            prediction_value = prediction_dict[threshold.name]
            prediction = assess_prediction_tool(threshold, prediction_value)
        except KeyError:
            prediction = None
        clinvar_same_nucleotide = clinvar_result[ClinVar_Type.SAME_NUCLEOTIDE]
        if clinvar_same_nucleotide.pathogenic and prediction:
            comment = f"The following ClinVar entries show splice variants at the same nucleotide position to be pathogenic: {clinvar_same_nucleotide.ids}."
            result = True
        elif not prediction:
            comment = f"The variant is not predicted to affect splicing."
            result = False
        else:
            comment = "No ClinVar entries found that show splice variants at the same nucleotide position as pathogenic.."
            result = False
        return RuleResult(
            "PS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )
