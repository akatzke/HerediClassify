#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    abstract_rule,
    rule_type,
)
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from information import Classification_Info, Info
from var_type import VARTYPE_GROUPS
from variant import TranscriptInfo


class Pp3_protein(abstract_rule):
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
                class_info.THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC,
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
        elif num_thresholds_met > 0:
            comment = f"Variant is predicted to be pathogenic by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
            result = True
        else:
            comment = f"Variant is not predicted to be pathogenic by {threshold.name} (threshold: {threshold.thresholds[0]}, value: {prediction_value})."
            result = False
        return RuleResult(
            "PP3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Pp3_splicing(abstract_rule):
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
        elif num_thresholds_met > 0:
            comment = f"Variant is predicted to have a splice effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
            result = True
        else:
            comment = f"Variant is not predicted to have a splice effect by {threshold.name} (threshold: {threshold.thresholds[0]}, value: {prediction_value})."
            result = False
        return RuleResult(
            "PP3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Pp3_splicing_cdh1(abstract_rule):
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
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        transcripts: list[TranscriptInfo],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        variant_types = [
            var_type for transcript in transcripts for var_type in transcript.var_type
        ]
        if num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}."
            result = False
        elif any(
            var_type in VARTYPE_GROUPS.INTRONIC.value for var_type in variant_types
        ):
            comment = f"Variant is located within cannonical splice site. PP3 does not apply here."
            result = False
        elif num_thresholds_met > 0:
            comment = f"Variant is predicted to have a splice effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
            result = True
        else:
            comment = f"Variant is not predicted to have a splice effect by {threshold.name} (threshold: {threshold.thresholds[0]}, value: {prediction_value})."
            result = False
        return RuleResult(
            "PP3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
