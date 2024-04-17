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
from acmg_rules.computation_evidence_utils import (
    assess_thresholds,
    Threshold,
)
from variant import TranscriptInfo, VariantInfo


class Bp4_protein(abstract_rule):
    """
    BP4: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_PATHOGENICITY_PREDICTION_BENIGN,
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
            comment = f"Variant is predicted to be benign by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
            result = True
        else:
            comment = f"Variant is not predicted to be benign by {threshold.name} (threshold: {threshold.thresholds[0]}, value: {prediction_value})."
            result = False
        return RuleResult(
            "BP4",
            rule_type.PROTEIN,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Bp4_protein_enigma(abstract_rule):
    """
    BP4: Assess results of prediction programs
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
                class_info.THRESHOLD_PATHOGENICITY_PREDICTION_BENIGN,
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
        if not any(
            type in VARTYPE_GROUPS.PREDICTION_PROTEIN.value for type in variant_types
        ):
            result = False
            comment = f"BP4 does not apply to this variant, as BP4 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
        elif num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}"
            result = False
        elif num_thresholds_met > 0:
            comment = f"Variant is predicted to be benign by {threshold.name}."
            result = True
        else:
            comment = f"Variant is not predicted to be benign by {threshold.name}."
            result = False
        return RuleResult(
            "BP4",
            rule_type.PROTEIN,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Bp4_splicing(abstract_rule):
    """
    BP4: Assess results of prediction programs
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
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
            comment = f"Variant is predicted to have no splice effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
            result = True
        else:
            comment = f"Variant is not predicted to have no splice effect by {threshold.name} (threshold: {threshold.thresholds[0]}, value: {prediction_value})."
            result = False
        return RuleResult(
            "BP4",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Bp4_splicing_enigma(abstract_rule):
    """
    BP4: Assess results of prediction programs
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
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
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
        if not any(
            type in VARTYPE_GROUPS.PREDICTION_SPLICING.value for type in variant_types
        ) and not any(
            type in VARTYPE_GROUPS.PREDICTION_PROTEIN.value for type in variant_types
        ):
            result = False
            comment = f"BP4 does not apply to this variant, as BP4 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
        elif num_thresholds_met is None:
            comment = f"No score was provided for {threshold.name}"
            result = False
        elif num_thresholds_met > 0:
            comment = (
                f"Variant is predicted to have no splicing effect by {threshold.name}."
            )
            result = True
        else:
            comment = f"Variant is not predicted to have no splicing effect by {threshold.name}."
            result = False
        return RuleResult(
            "BP4",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
