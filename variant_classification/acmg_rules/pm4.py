#!/usr/bin/env python3

from typing import Callable

import pathlib

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    abstract_rule,
    rule_type,
    summarise_results_per_transcript,
)
from information import Info, Classification_Info
from transcript_annotated import TranscriptInfo_annot
from variant import VariantInfo
from var_type import VARTYPE
from utils import select_mane_transcript


class Pm4(abstract_rule):
    """
    PM4: Protein length change caused by variant is above 10% threshold
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
                class_info.MANE_TRANSCRIPT_LIST_PATH,
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript_list: list[TranscriptInfo_annot],
        variant: VariantInfo,
        mane_path: pathlib.Path,
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        if not annotated_transcript_list:
            return RuleResult(
                "PM4",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                f"No annotated transcripts provided, PM4 can not be applied",
            )
        if len(annotated_transcript_list):
            transcript = annotated_transcript_list[0]
        else:
            transcript = select_mane_transcript(annotated_transcript_list, mane_path)
        if transcript is None:
            comment = f"PM4 does not apply to this variant, as PM4 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            result = False
        elif (
            transcript.diff_len_protein_percent > threshold_diff_len_prot_percent
            and not transcript.len_change_in_repetitive_region
        ):
            comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. A repetitive region is not affected."
            result = True
        elif (
            transcript.diff_len_protein_percent > threshold_diff_len_prot_percent
            and transcript.len_change_in_repetitive_region
        ):
            comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region overlaps repetitive region."
            result = False
        else:
            comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}."
            result = False
        return RuleResult(
            "PM4",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )


class Pm4_stoploss(abstract_rule):
    """
    PM4: Protein length change caused by variant is above 10% threshold
    Only applicable to stop loss variants
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.ANNOTATED_TRANSCRIPT_LIST,
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
                class_info.VARIANT,
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript_list: list[TranscriptInfo_annot],
        threshold_diff_len_prot_percent: float,
        variant: VariantInfo,
        mane_path: pathlib.Path,
    ) -> RuleResult:
        results = {}
        for transcript in annotated_transcript_list:
            if any(var_type is VARTYPE.STOP_LOST for var_type in transcript.var_type):
                if (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                    and transcript.is_truncated_region_disease_relevant
                    and not transcript.len_change_in_repetitive_region
                ):
                    comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region does not overlap repetitive region."
                    result = True
                elif (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                    and transcript.is_truncated_region_disease_relevant
                    and transcript.len_change_in_repetitive_region
                ):
                    comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region overlaps repetitive region."
                    result = False
                else:
                    comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}"
                    result = False
                rule_result = RuleResult(
                    "PM4",
                    rule_type.GENERAL,
                    evidence_type.PATHOGENIC,
                    result,
                    evidence_strength.MODERATE,
                    comment,
                )
                results[transcript.transcript_id] = rule_result
        if len(results) == 0:
            comment = f"PM4 does not apply to this variant, as PM4 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            final_result = RuleResult(
                "PM4",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                comment,
            )
        else:
            final_result = summarise_results_per_transcript(results, mane_path)
        return final_result


class Pm4_pten(abstract_rule):
    """
    PM4: Protein length change caused by variant is above 10% threshold
    Only applicable to stop loss variants
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
                class_info.VARIANT_HOTSPOT_ANNOTATION,
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript_list: list[TranscriptInfo_annot],
        variant: VariantInfo,
        variant_in_hotspot: bool,
        mane_path: pathlib.Path,
    ) -> RuleResult:
        results = {}
        for transcript in annotated_transcript_list:
            # All extensions of the protein can be seen as pathogenic
            if any(var_type is VARTYPE.STOP_LOST for var_type in transcript.var_type):
                if not transcript.len_change_in_repetitive_region:
                    comment = f"Length of disease relevant transcript {transcript.transcript_id} is increased by {abs(transcript.diff_len_protein_percent)}. Deleted region does not overlap repetitive region."
                    result = True
                else:
                    comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}"
                    result = False
                rule_result = RuleResult(
                    "PM4",
                    rule_type.GENERAL,
                    evidence_type.PATHOGENIC,
                    result,
                    evidence_strength.MODERATE,
                    comment,
                )
                results[transcript.transcript_id] = rule_result
            elif any(
                var_type in [VARTYPE.INFRAME_DELETION, VARTYPE.INFRAME_INSERTION]
                for var_type in transcript.var_type
            ):
                if (
                    not transcript.len_change_in_repetitive_region
                    and variant_in_hotspot
                ):
                    comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region is overlaps mutational hotspot."
                    result = True
                elif transcript.len_change_in_repetitive_region:
                    comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region overlaps repetitive region."
                    result = False
                else:
                    comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}"
                    result = False
                rule_result = RuleResult(
                    "PM4",
                    rule_type.GENERAL,
                    evidence_type.PATHOGENIC,
                    result,
                    evidence_strength.MODERATE,
                    comment,
                )
                results[transcript.transcript_id] = rule_result
            elif any(
                var_type is VARTYPE.FRAMESHIFT_VARIANT
                for var_type in transcript.var_type
            ):
                if transcript.diff_len_protein_percent > 0:
                    result = False
                    comment = f"Frameshift variant decreases the length of the transcript, PM4 only applies to extensions of the protein."
                elif not transcript.len_change_in_repetitive_region:
                    comment = f"Length of disease relevant transcript {transcript.transcript_id} is increased by {abs(transcript.diff_len_protein_percent)}."
                    result = True
                else:
                    comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}"
                    result = False
                rule_result = RuleResult(
                    "PM4",
                    rule_type.GENERAL,
                    evidence_type.PATHOGENIC,
                    result,
                    evidence_strength.MODERATE,
                    comment,
                )
                results[transcript.transcript_id] = rule_result

        if len(results) == 0:
            comment = f"PM4 does not apply to this variant, as PM4 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            final_result = RuleResult(
                "PM4",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                comment,
            )
        else:
            final_result = summarise_results_per_transcript(results, mane_path)
        return final_result
