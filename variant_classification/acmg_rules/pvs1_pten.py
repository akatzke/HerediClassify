#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)
from acmg_rules.pvs1 import Pvs1
from information import Classification_Info, Info
from variant import RNAData, TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)
from acmg_rules.computation_evidence_utils import (
    Threshold,
    assess_thresholds,
)
from acmg_rules.functional_splicing_assay_utils import (
    adjust_strength_according_to_rna_data_pvs1,
)


class Pvs1_pten(Pvs1):
    """
    PVS1: loss of function
    Following VCEP guidelines for PTEN
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.ANNOTATED_TRANSCRIPT_LIST,
                class_info.VARIANT,
                class_info.SPLICING_ASSAY,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcripts: list[TranscriptInfo],
        variant: VariantInfo,
        splice_assay: Optional[list[RNAData]],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ):
        results = []
        for transcript in annotated_transcripts:
            if isinstance(transcript, TranscriptInfo_exonic):
                result = cls.assess_pvs1_frameshift_PTC_pten(transcript)
                results.append(result)
            elif isinstance(transcript, TranscriptInfo_intronic):
                splice_result = cls.assess_pvs1_splice_pten(
                    transcript, prediction_dict, threshold
                )
                if splice_assay:
                    splice_result = adjust_strength_according_to_rna_data_pvs1(
                        splice_assay, splice_result
                    )
                results.append(splice_result)
            elif isinstance(transcript, TranscriptInfo_start_loss):
                result = cls.assess_pvs1_start_loss_pathogenic_very_strong()
                results.append(result)
        if len(results) == 0:
            if not annotated_transcripts:
                comment = "No annotated transcripts provided, PVS1 can not be applied."
            else:
                comment = f"PVS1 does not apply to this variant, as PVS1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            result = RuleResult(
                "PVS1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                comment,
            )
        else:
            result = summarise_results_per_transcript(results)
        return result

    @classmethod
    def assess_pvs1_frameshift_PTC_pten(
        cls, transcript: TranscriptInfo_exonic
    ) -> RuleResult:
        if transcript.is_NMD and not transcript.is_reading_frame_preserved:
            comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD and in a disease relevant transcript."
            result = True
            strength = evidence_strength.VERY_STRONG
        elif (
            transcript.is_truncated_region_disease_relevant
            and transcript.is_reading_frame_preserved
        ):
            comment = (
                f"Transcript {transcript.transcript_id} is not predict to undergo NMD. Truncated region is disease relevant. "
                + transcript.comment_truncated_region
            )
            result = True
            strength = evidence_strength.MODERATE
        else:
            comment = f"Variant in transcript {transcript.transcript_id} does not meet any of the specified criteria for PTEN."
            result = False
            strength = evidence_strength.VERY_STRONG
        return RuleResult(
            "PVS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )

    @classmethod
    def assess_pvs1_splice_pten(
        cls,
        transcript: TranscriptInfo_intronic,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if not transcript.are_exons_skipped and not num_thresholds_met:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"No splicing alteration predicted for transcript {transcript.transcript_id}."
        elif not transcript.is_reading_frame_preserved and transcript.is_NMD:
            comment = f"Transcript {transcript.transcript_id} is predicted undergo NMD and is disease relevant"
            result = True
            strength = evidence_strength.VERY_STRONG
        elif (
            transcript.is_reading_frame_preserved
            and transcript.is_truncated_region_disease_relevant
        ):
            comment = (
                f"Transcript {transcript.transcript_id}'s reading frame is preserved and skipped exon is disease relevant. "
                + transcript.comment_truncated_region
            )
            result = True
            strength = evidence_strength.STRONG
        else:
            comment = f"Variant in transcript {transcript.transcript_id} does not meet any of the specified criteria for PTEN."
            result = False
            strength = evidence_strength.VERY_STRONG
        return RuleResult(
            "PVS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
