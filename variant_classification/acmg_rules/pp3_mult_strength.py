#!/usr/bin/env python3

from typing import Callable

import pathlib

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    abstract_rule,
    rule_type,
)
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from var_type import VARTYPE_GROUPS
from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo
from utils import check_intersection_with_bed
from ensembl import ensembl


class Pp3_protein_mult_strength(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.TRANSCRIPT,
                class_info.VARIANT,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        transcripts: list[TranscriptInfo],
        variant: VariantInfo,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        if len(transcripts) == 1:
            variant_types = transcripts[0].var_type
        else:
            variant_types = variant.var_type
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if any(
            type in VARTYPE_GROUPS.PREDICTION_NO_PROTEIN.value for type in variant_types
        ):
            comment = f"PP3 does not apply to this variant, as PP3 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
            strength = evidence_strength.SUPPORTING
            result = False
        elif num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}."
            strength = evidence_strength.SUPPORTING
            result = False
        elif num_thresholds_met == 0:
            comment = f"Variant is not predicted to be pathogenic by {threshold.name}."
            strength = evidence_strength.SUPPORTING
            result = False
        else:
            result = True
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to be pathogenic by {threshold.name} with evidence strength {strength.value} meeting a threshold of {threshold.thresholds[num_thresholds_met -1]} (value: {prediction_value})."
        return RuleResult(
            "PP3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )


class Pp3_splicing_mult_strength(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}."
            result = False
            strength = evidence_strength.SUPPORTING
        elif num_thresholds_met == 0:
            comment = (
                f"Variant is not predicted to have a splice effect by {threshold.name}."
            )
            result = False
            strength = evidence_strength.SUPPORTING
        else:
            result = True
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to have a splice effect by {threshold.name} with evidence strength {strength.value} meeting a threshold of {threshold.thresholds[num_thresholds_met -1]} (value: {prediction_value})."
        return RuleResult(
            "PP3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )


class Pp3_splicing_enigma_mult_strength(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.TRANSCRIPT,
                class_info.VARIANT,
                class_info.CRITICAL_REGION_PATH,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        transcripts: list[TranscriptInfo],
        variant: VariantInfo,
        path_critical_region: pathlib.Path,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        if len(transcripts) == 1:
            variant_types = transcripts[0].var_type
        else:
            variant_types = variant.var_type
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if not transcripts:
            return RuleResult(
                "PP3",
                rule_type.SPLICING,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"Variant is not coding in disease relevant transcript and is therfore located outside of the disease relevant region defined by the VCEP.",
            )
        ref_transcript = ensembl.transcript_by_id(transcripts[0].transcript_id)
        in_critical_region, _ = check_intersection_with_bed(
            variant,
            variant.genomic_start,
            variant.genomic_end,
            ref_transcript,
            path_critical_region,
        )
        if not any(
            type in VARTYPE_GROUPS.PREDICTION_SPLICING.value for type in variant_types
        ):
            strength = evidence_strength.SUPPORTING
            result = False
            comment = f"PP3 does not apply to this variant, as PP3 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
        elif not in_critical_region:
            strength = evidence_strength.SUPPORTING
            result = False
            comment = f"PP3 does not apply to this variant, as the variant is located outside of the disease relevant region defined by the VCEP."
        elif num_thresholds_met is None:
            strength = evidence_strength.SUPPORTING
            comment = f"No score was provided for {threshold.name}."
            result = False
        elif num_thresholds_met > 0:
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to have a splice effect by {threshold.name} with evidence strength {strength.value} meeting a threshold of {threshold.thresholds[num_thresholds_met -1]} (value: {prediction_value})."
            result = True
        else:
            strength = evidence_strength.SUPPORTING
            comment = f"Variant is not predicted to have a splice effect by {threshold.name} (threshold: {threshold.thresholds[0]}, value: {prediction_value})."
            result = False
        return RuleResult(
            "PP3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )


class Pp3_protein_enigma_mult_strength(abstract_rule):
    """
    PP3: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.TRANSCRIPT,
                class_info.VARIANT,
                class_info.CRITICAL_REGION_PATH,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        transcripts: list[TranscriptInfo],
        variant: VariantInfo,
        path_critical_region: pathlib.Path,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        if len(transcripts) == 1:
            variant_types = transcripts[0].var_type
        else:
            variant_types = variant.var_type
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if not transcripts:
            return RuleResult(
                "PP3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"Variant is not coding in disease relevant transcript and is therfore located outside of the disease relevant region defined by the VCEP.",
            )
        ref_transcript = ensembl.transcript_by_id(transcripts[0].transcript_id)
        in_critical_region, _ = check_intersection_with_bed(
            variant,
            variant.genomic_start,
            variant.genomic_end,
            ref_transcript,
            path_critical_region,
        )
        if not any(
            type in VARTYPE_GROUPS.PREDICTION_PROTEIN.value for type in variant_types
        ):
            strength = evidence_strength.SUPPORTING
            result = False
            comment = f"PP3 does not apply to this variant, as PP3 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
        elif not in_critical_region:
            strength = evidence_strength.SUPPORTING
            result = False
            comment = f"PP3 does not apply to this variant, as the variant is located outside of the disease relevant region defined by the VCEP."
        elif num_thresholds_met is None:
            strength = evidence_strength.SUPPORTING
            comment = f"No score was provided for {threshold.name}."
            result = False
        elif num_thresholds_met > 0:
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Variant is predicted to be pathogenic by {threshold.name} with evidence strength {strength.value} meeting a threshold of {threshold.thresholds[num_thresholds_met -1]} (value: {prediction_value})."
            result = True
        else:
            strength = evidence_strength.SUPPORTING
            comment = f"Variant is not predicted to be pathogenic by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met -1]}, value: {prediction_value})."
            result = False
        return RuleResult(
            "PP3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
