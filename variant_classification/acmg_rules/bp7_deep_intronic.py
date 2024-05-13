#!/usr/bin/env python3

import pathlib

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info
from acmg_rules.computation_evidence_utils import (
    assess_thresholds,
    Threshold,
)
from acmg_rules.functional_splicing_assay_utils import (
    assess_splicing_data_bp7,
)
from variant import RNAData, TranscriptInfo, VariantInfo
from var_type import VARTYPE
from utils import select_mane_transcript


class Bp7_deep_intronic_enigma(abstract_rule):
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
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        transcripts: list[TranscriptInfo],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splicing_assay: Optional[list[RNAData]],
        mane_path: pathlib.Path,
    ) -> RuleResult:
        # Check RNA data first
        if splicing_assay:
            performed, result_assay, comment_assay = assess_splicing_data_bp7(
                splicing_assay
            )
            if performed:
                return RuleResult(
                    "BP7_RNA",
                    rule_type.SPLICING,
                    evidence_type.BENIGN,
                    result_assay,
                    evidence_strength.STRONG,
                    comment_assay,
                )

        # In case one disase variant transcripts is defined, use type of variant in that transcript
        # Otherwise use all variant types defined for variant
        if len(transcripts) == 1:
            transcript = transcripts[0]
            variant_types = transcript.var_type
        else:
            transcript = select_mane_transcript(transcripts, mane_path)
            if transcript is None:
                variant_types = variant.var_type
            else:
                variant_types = transcript.var_type
        # Check prediction
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
            result = False
            comment = f"No score was provided for {threshold.name}"
        elif num_thresholds_met > 0:
            """
            Variant is predicted to have splicing effect
            Now check if variant type applies
            """
            if any(
                var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant_types
            ):
                result = True
                comment = f"The synonymous variant is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                if transcript is None:
                    result = False
                    comment = (
                        comment
                        + " The variant is not located in the MANE transcript. Therefore BP7 is disabled."
                    )

            elif any(var_type is VARTYPE.INTRON_VARIANT for var_type in variant_types):
                """
                Check if variant is a deep intronic variant
                """
                result = False
                comments_all = []
                for trans in transcripts:
                    if (
                        trans.var_hgvs.pos.start.offset >= 7
                        and trans.var_hgvs.pos.end.offset >= 7
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {trans.transcript_id} is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                        comments_all.append(comment_tmp)
                    elif (
                        trans.var_hgvs.pos.start.offset <= -21
                        and trans.var_hgvs.pos.end.offset <= -21
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {trans.transcript_id} is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                        comments_all.append(comment_tmp)
                if result:
                    comment = " ".join(comments_all)
                else:
                    comment = f"BP7 does not apply to this variant, as it is not located in the defined region for deep intronic region variants (<= -21 or >= 7)."
                if transcript is None:
                    result = False
                    comment = (
                        comment
                        + " The variant is not located in the MANE transcript. Therefore BP7 is disabled."
                    )
            else:
                result = False
                comment = f"BP7 does not apply to this variant, as BP7 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
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


class Bp7_deep_intronic_atm(abstract_rule):
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
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        transcripts: list[TranscriptInfo],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splicing_assay: Optional[list[RNAData]],
        mane_path: pathlib.Path,
    ) -> RuleResult:
        # Check RNA data first
        if splicing_assay:
            performed, result_assay, comment_assay = assess_splicing_data_bp7(
                splicing_assay
            )
            if performed:
                return RuleResult(
                    "BP7_RNA",
                    rule_type.SPLICING,
                    evidence_type.BENIGN,
                    result_assay,
                    evidence_strength.STRONG,
                    comment_assay,
                )

        # In case one disase variant transcripts is defined, use type of variant in that transcript
        # Otherwise use all variant types defined for variant
        if len(transcripts) == 1:
            transcript = transcripts[0]
            variant_types = transcript.var_type
        else:
            transcript = select_mane_transcript(transcripts, mane_path)
            if transcript is None:
                variant_types = variant.var_type
            else:
                variant_types = transcript.var_type
        # Check prediction
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
            result = False
            comment = f"No score was provided for {threshold.name}"
        elif num_thresholds_met > 0:
            """
            Variant is predicted to have splicing effect
            Now check if variant type applies
            """
            if any(
                var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant_types
            ):
                result = True
                comment = f"The synonymous variant is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                if transcript is None:
                    result = False
                    comment = (
                        comment
                        + " The variant is not located in the MANE transcript. Therefore BP7 is disabled."
                    )

            elif any(var_type is VARTYPE.INTRON_VARIANT for var_type in variant_types):
                """
                Check if variant is a deep intronic variant
                """
                result = False
                comments_all = []
                for trans in transcripts:
                    if (
                        trans.var_hgvs.pos.start.offset > 7
                        and trans.var_hgvs.pos.end.offset > 7
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {trans.transcript_id} is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                        comments_all.append(comment_tmp)
                    elif (
                        trans.var_hgvs.pos.start.offset < -40
                        and trans.var_hgvs.pos.end.offset < -40
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {trans.transcript_id} is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                        comments_all.append(comment_tmp)
                if result:
                    comment = " ".join(comments_all)
                else:
                    comment = f"BP7 does not apply to this variant, as it is not located in the defined region for deep intronic region variants (< -40 or > 7)."
                if transcript is None:
                    result = False
                    comment = (
                        comment
                        + " The variant is not located in the MANE transcript. Therefore BP7 is disabled."
                    )
            else:
                result = False
                comment = f"BP7 does not apply to this variant, as BP7 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
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


class Bp7_deep_intronic_palb2(abstract_rule):
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
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        transcripts: list[TranscriptInfo],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splicing_assay: Optional[list[RNAData]],
        mane_path: pathlib.Path,
    ) -> RuleResult:
        # Check RNA data first
        if splicing_assay:
            performed, result_assay, comment_assay = assess_splicing_data_bp7(
                splicing_assay
            )
            if performed:
                return RuleResult(
                    "BP7_RNA",
                    rule_type.SPLICING,
                    evidence_type.BENIGN,
                    result_assay,
                    evidence_strength.STRONG,
                    comment_assay,
                )

        # In case one disase variant transcripts is defined, use type of variant in that transcript
        # Otherwise use all variant types defined for variant
        if len(transcripts) == 1:
            transcript = transcripts[0]
            variant_types = transcript.var_type
        else:
            transcript = select_mane_transcript(transcripts, mane_path)
            if transcript is None:
                variant_types = variant.var_type
            else:
                variant_types = transcript.var_type
        # Check prediction
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
            result = False
            comment = f"No score was provided for {threshold.name}"
        elif num_thresholds_met > 0:
            """
            Variant is predicted to have splicing effect
            Now check if variant type applies
            """
            if any(
                var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant_types
            ):
                result = True
                comment = f"The synonymous variant is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                if transcript is None:
                    result = False
                    comment = (
                        comment
                        + " The variant is not located in the MANE transcript. Therefore BP7 is disabled."
                    )

            elif any(var_type is VARTYPE.INTRON_VARIANT for var_type in variant_types):
                """
                Check if variant is a deep intronic variant
                """
                result = False
                comments_all = []
                for trans in transcripts:
                    if (
                        trans.var_hgvs.pos.start.offset > 7
                        and trans.var_hgvs.pos.end.offset > 7
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {trans.transcript_id} is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                        comments_all.append(comment_tmp)
                    elif (
                        trans.var_hgvs.pos.start.offset < -21
                        and trans.var_hgvs.pos.end.offset < -21
                    ):
                        result = True
                        comment_tmp = f"The deep intronic variant in {trans.transcript_id} is predicted to have no splicing effect by {threshold.name} (threshold: {threshold.thresholds[num_thresholds_met-1]}, value: {prediction_value})."
                        comments_all.append(comment_tmp)
                if result:
                    comment = " ".join(comments_all)
                else:
                    comment = f"BP7 does not apply to this variant, as it is not located in the defined region for deep intronic region variants (< -21 or > 7)."
                if transcript is None:
                    result = False
                    comment = (
                        comment
                        + " The variant is not located in the MANE transcript. Therefore BP7 is disabled."
                    )
            else:
                result = False
                comment = f"BP7 does not apply to this variant, as BP7 does not apply to variant types {', '.join([var_type.value for var_type in variant_types])}."
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
