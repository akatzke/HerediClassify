#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info
from acmg_rules.computation_evidence_utils import (
    assess_prediction_tool,
    Threshold,
)
from variant import FunctionalData, TranscriptInfo, VariantInfo
from var_type import VARTYPE


class Bp7_with_rna_assay_deep_intronic_enigma(abstract_rule):
    """
    BP7: Silent missense variant is predicted to have effect on splicing
    Expanded to include deep intronic variants
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.TRANSCRIPT,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
                class_info.SPLICING_ASSAY,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        transcripts: list[TranscriptInfo],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splicing_assay: FunctionalData,
    ) -> RuleResult:
        if splicing_assay.performed:
            if splicing_assay.benign:
                comment = f"Splicing assay shows no effect on splicing."
                result = True
            else:
                comment = f"Splicing assay does not show no effect on splicing."
                result = False
            return RuleResult(
                "BP7",
                rule_type.SPLICING,
                evidence_type.BENIGN,
                result,
                evidence_strength.STRONG,
                comment,
            )

        prediction_value = prediction_dict.get(threshold.name, None)
        prediction = assess_prediction_tool(threshold, prediction_value)
        if prediction is None:
            result = False
            comment = f"No score was provided for {threshold.name}"
        elif prediction:
            """
            Variant is predicted to have splicing effect
            Now check if variant type applies
            """
            if any(
                var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant.var_type
            ):
                result = True
                comment = f"The synonymous variant is predicted to have no splicing effect by {threshold.name}."

            else:
                """
                Check if variant is a deep intronic variant
                """
                prediction_value = prediction_dict.get(threshold.name, None)
                prediction = assess_prediction_tool(threshold, prediction_value)
                result = False
                comments_all = []
                for transcript in transcripts:
                    if (
                        transcript.var_hgvs.pos.start.offset >= 7
                        and transcript.var_hgvs.pos.end.offset >= 7
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {transcript.transcript_id} is predicted to have no splicing effect by {threshold.name}."
                        comments_all.append(comment_tmp)
                    elif (
                        transcript.var_hgvs.pos.start.offset <= -21
                        and transcript.var_hgvs.pos.end.offset <= -21
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {transcript.transcript_id} is predicted to have no splicing effect by {threshold.name}."
                        comments_all.append(comment_tmp)
                if result:
                    comment = " ".join(comments_all)
                else:
                    comment = f"BP7 does not apply to this variant, as BP7 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
        else:
            result = False
            comment = f"The variant is not predicted to not affect splicing by {threshold.name}."

        return RuleResult(
            "BP7",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Bp7_with_rna_assay_deep_intronic_atm(abstract_rule):
    """
    BP7: Silent missense variant is predicted to have effect on splicing
    Expanded to include deep intronic variants
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.TRANSCRIPT,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
                class_info.SPLICING_ASSAY,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        transcripts: list[TranscriptInfo],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splicing_assay: FunctionalData,
    ) -> RuleResult:
        if splicing_assay.performed:
            if splicing_assay.benign:
                comment = f"Splicing assay shows no effect on splicing."
                result = True
            else:
                comment = f"Splicing assay does not show no effect on splicing."
                result = False
            return RuleResult(
                "BP7",
                rule_type.SPLICING,
                evidence_type.BENIGN,
                result,
                evidence_strength.STRONG,
                comment,
            )

        prediction_value = prediction_dict.get(threshold.name, None)
        prediction = assess_prediction_tool(threshold, prediction_value)
        if prediction is None:
            result = False
            comment = f"No score was provided for {threshold.name}"
        elif prediction:
            """
            Variant is predicted to have splicing effect
            Now check if variant type applies
            """
            if any(
                var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant.var_type
            ):
                result = True
                comment = f"The synonymous variant is predicted to have no splicing effect by {threshold.name}."

            else:
                """
                Check if variant is a deep intronic variant
                """
                prediction_value = prediction_dict.get(threshold.name, None)
                prediction = assess_prediction_tool(threshold, prediction_value)
                result = False
                comments_all = []
                for transcript in transcripts:
                    if (
                        transcript.var_hgvs.pos.start.offset > 7
                        and transcript.var_hgvs.pos.end.offset > 7
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {transcript.transcript_id} is predicted to have no splicing effect by {threshold.name}."
                        comments_all.append(comment_tmp)
                    elif (
                        transcript.var_hgvs.pos.start.offset < -40
                        and transcript.var_hgvs.pos.end.offset < -40
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {transcript.transcript_id} is predicted to have no splicing effect by {threshold.name}."
                        comments_all.append(comment_tmp)
                if result:
                    comment = " ".join(comments_all)
                else:
                    comment = f"BP7 does not apply to this variant, as BP7 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
        else:
            result = False
            comment = f"The variant is not predicted to not affect splicing by {threshold.name}."

        return RuleResult(
            "BP7",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Bp7_with_rna_assay_deep_intronic_palb2(abstract_rule):
    """
    BP7: Silent missense variant is predicted to have effect on splicing
    Expanded to include deep intronic variants
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.TRANSCRIPT,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
                class_info.SPLICING_ASSAY,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        transcripts: list[TranscriptInfo],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splicing_assay: FunctionalData,
    ) -> RuleResult:
        if splicing_assay.performed:
            if splicing_assay.benign:
                comment = f"Splicing assay shows no effect on splicing."
                result = True
            else:
                comment = f"Splicing assay does not show no effect on splicing."
                result = False
            return RuleResult(
                "BP7",
                rule_type.SPLICING,
                evidence_type.BENIGN,
                result,
                evidence_strength.STRONG,
                comment,
            )

        prediction_value = prediction_dict.get(threshold.name, None)
        prediction = assess_prediction_tool(threshold, prediction_value)
        if prediction is None:
            result = False
            comment = f"No score was provided for {threshold.name}"
        elif prediction:
            """
            Variant is predicted to have splicing effect
            Now check if variant type applies
            """
            if any(
                var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant.var_type
            ):
                result = True
                comment = f"The synonymous variant is predicted to have no splicing effect by {threshold.name}."

            else:
                """
                Check if variant is a deep intronic variant
                """
                prediction_value = prediction_dict.get(threshold.name, None)
                prediction = assess_prediction_tool(threshold, prediction_value)
                result = False
                comments_all = []
                for transcript in transcripts:
                    if (
                        transcript.var_hgvs.pos.start.offset > 7
                        and transcript.var_hgvs.pos.end.offset > 7
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {transcript.transcript_id} is predicted to have no splicing effect by {threshold.name}."
                        comments_all.append(comment_tmp)
                    elif (
                        transcript.var_hgvs.pos.start.offset < -21
                        and transcript.var_hgvs.pos.end.offset < -21
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {transcript.transcript_id} is predicted to have no splicing effect by {threshold.name}."
                        comments_all.append(comment_tmp)
                if result:
                    comment = " ".join(comments_all)
                else:
                    comment = f"BP7 does not apply to this variant, as BP7 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
        else:
            result = False
            comment = f"The variant is not predicted to not affect splicing by {threshold.name}."

        return RuleResult(
            "BP7",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
