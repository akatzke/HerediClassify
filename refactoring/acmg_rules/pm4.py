#!/usr/bin/env python3
from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    abstract_rule,
    rule_type,
    summarise_results_per_transcript,
)
from refactoring.information import Info, Classification_Info
from refactoring.transcript_annotated import TranscriptInfo_annot


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
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript_list: list[TranscriptInfo_annot],
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        results = []
        for transcript in annotated_transcript_list:
            if not transcript.is_truncated_region_disease_relevant:
                comment = (
                    f"Transcript {transcript.transcript_id} is not disease relevant."
                )
                result = False
            elif (
                transcript.diff_len_protein_percent > threshold_diff_len_prot_percent
                and transcript.is_truncated_region_disease_relevant
                and not transcript.len_change_in_repetitive_region
            ):
                comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region does not overlap repetitive region."
                result = True
            elif (
                transcript.diff_len_protein_percent > 0.1
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
            results.append(rule_result)
        final_result = summarise_results_per_transcript(results)
        return final_result
