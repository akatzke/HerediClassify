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
from variant import RNAData, TranscriptInfo, VariantInfo
from var_type import VARTYPE
from acmg_rules.functional_splicing_assay_utils import (
    assess_splicing_data_bp7,
)
from utils import select_mane_transcript


class Bp7(abstract_rule):
    """
    BP7: Silent missense variant is predicted to have effect on splicing
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
                class_info.SPLICING_ASSAY,
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        transcripts: list[TranscriptInfo],
        variant: VariantInfo,
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
            variant_types = transcripts[0].var_type
        else:
            transcript = select_mane_transcript(transcripts, mane_path)
            variant_types = transcript.var_type

        # Check prediction
        if not any(
            var_type is VARTYPE.SYNONYMOUS_VARIANT for var_type in variant_types
        ):
            return RuleResult(
                "BP7",
                rule_type.SPLICING,
                evidence_type.BENIGN,
                False,
                evidence_strength.SUPPORTING,
                f"BP7 does not apply to this variant, as BP7 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}.",
            )
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if num_thresholds_met is None:
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
            "BP7",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )
