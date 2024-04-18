#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from var_type import VARTYPE_GROUPS
from information import Classification_Info, Info
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from clinvar_utils import ClinVar_Status, ClinVar_Type, ClinVar, get_affected_transcript
from variant import TranscriptInfo, VariantInfo
from custom_exceptions import No_transcript_with_var_type_found


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
            if clinvar_same_aa.associated_ids:
                comment = (
                    comment
                    + f" The following ClinVar entries show  same amino acid change as likely pathogenic: {clinvar_same_aa.associated_ids}."
                )
        else:
            comment = "No ClinVar entries found that show the same amino acid change as pathogenic."
            result = False
        return RuleResult(
            "PS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Ps1_protein_spliceai(abstract_rule):
    """
    PS1: Position is classified as pathogenic
    The splicing effect of the missense variant is excluded based on the SpliceAI prediction
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
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_result: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        clinvar_same_aa = clinvar_result[ClinVar_Type.SAME_AA_CHANGE]
        if clinvar_same_aa.pathogenic:
            if num_thresholds_met is None:
                result = True
                comment = f"ATTENTION: No splicing prediction is available for variant under assessment. The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}."
            elif num_thresholds_met > 0:
                comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}."
                result = True
                if clinvar_same_aa.associated_ids:
                    comment = (
                        comment
                        + f" The following ClinVar entreis show the same amino acid exchange as likely pathogenic: {clinvar_same_aa.associated_ids}."
                    )
            else:
                comment = f"Variant is not predicted to not affect splicing. PS1_protein is therehfore not applicable."
                result = False
        else:
            comment = "No ClinVar entries found that show the same amino acid change as pathogneic."
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
            if clinvar_same_nucleotide.associated_ids:
                comment = (
                    comment
                    + f" The following ClinVar entreis show splice variants at the same nucleotide position as likely pathogenic: {clinvar_same_nucleotide.associated_ids}."
                )
        else:
            comment = "No ClinVar entries found that show splice variants at the same nucleotide position as pathogenic."
            result = False
        return RuleResult(
            "PS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Ps1_splicing_clingen(abstract_rule):
    """
    PS1 for splicing: Splice variant in same position has been show to be pathogenic
    Following the recommendations of the ClinGen splicing group
    The adjustment of evidence strength in case PVS1 applies to the variant is handled in elsewhere
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_CLINVAR_SPLICEAI_SPLICE,
                class_info.TRANSCRIPT,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        clinvar_result: dict[ClinVar_Type, ClinVar],
        transcripts: list[TranscriptInfo],
    ) -> RuleResult:
        clinvar_same_nucleotide = clinvar_result[ClinVar_Type.SAME_NUCLEOTIDE]
        clinvar_same_splice_site = clinvar_result[ClinVar_Type.SAME_SPLICE_SITE]
        try:
            affected_transcript, _ = get_affected_transcript(
                transcripts, VARTYPE_GROUPS.INTRONIC
            )
        except No_transcript_with_var_type_found:
            return RuleResult(
                "PS1",
                rule_type.SPLICING,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.STRONG,
                f"PS1 splicing following ClinGen specifications does not apply to this variant, as PS1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}.",
            )

        if clinvar_same_nucleotide.pathogenic:
            result = True
            if (
                clinvar_same_nucleotide.highest_classification
                == ClinVar_Status.PATHOGENIC
            ):
                comment = f"The following ClinVar entries show splice variants at the same nucleotide position to be pathogenic: {clinvar_same_nucleotide.ids}."
                strength = evidence_strength.STRONG
                if clinvar_same_nucleotide.associated_ids:
                    comment = (
                        comment
                        + f" The following ClinVar entreis show splice variants at the same nucleotide position as likely pathogenic: {clinvar_same_nucleotide.associated_ids}."
                    )
            else:
                if (
                    abs(affected_transcript.var_hgvs.pos.start.offset) > 2
                    or abs(affected_transcript.var_hgvs.pos.start.offset) == 0
                ) and (
                    abs(affected_transcript.var_hgvs.pos.end.offset) > 2
                    or abs(affected_transcript.var_hgvs.pos.end.offset == 0)
                ):
                    # In case the canonical dinucleotide is not affected, PS1 supporting applies
                    comment = f"The following ClinVar entries show splice variants at the same nucleotide position to be likely pathogenic: {clinvar_same_nucleotide.ids}."
                    strength = evidence_strength.MODERATE
                else:
                    # In case the canonical dinucleotide is affected, this criterium does not apply and is set to false
                    result = False
                    strength = evidence_strength.STRONG
                    comment = f"The variant is located in the canonical dinucleotide and a likely pathogenic variant was found in ClinVar affecting the same splice site. PS1 is not applicable in this case."
        elif clinvar_same_splice_site.pathogenic:
            result = True
            if (
                clinvar_same_splice_site.highest_classification
                == ClinVar_Status.PATHOGENIC
            ):
                comment = f"The following ClinVar entries show splice variants at the same splice site to be pathogenic: {clinvar_same_splice_site.ids}."
                strength = evidence_strength.MODERATE
                if clinvar_same_splice_site.associated_ids:
                    comment = (
                        comment
                        + f" The following ClinVar entreis show splice variants at the same splice site as likely pathogenic: {clinvar_same_splice_site.associated_ids}."
                    )
            else:
                comment = f"The following ClinVar entries show splice variants at the same splice site to be likely pathogenic: {clinvar_same_splice_site.ids}."
                strength = evidence_strength.SUPPORTING
        else:
            result = False
            strength = evidence_strength.STRONG
            comment = "No ClinVar entries found that show splice variants at the same nucleotide position or at the same splice site as (likely) pathogenic."
        return RuleResult(
            "PS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
