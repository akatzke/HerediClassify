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
import refactoring.information as info
from refactoring.variant import TranscriptInfo
from refactoring.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
)


class bp3(abstract_rule):
    """
    BP3: Protein length change in repetitive region
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (
                info.classification_information.ANNOTATED_TRANSCRIPT_LIST,
                info.classification_information.THRESHOLD_DIFF_LEN_PROT_PERCENT,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcripts: list[TranscriptInfo],
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        results = []
        for transcript in annotated_transcripts:
            if (
                type(transcript) != TranscriptInfo_exonic
                or type(transcript) != TranscriptInfo_intronic
            ):
                comment = f"Transcript {transcript.transcript_id} does not carry variant of exonic or intronic variant type."
                result = False
                results.append(result)
                break
            if not transcript.transcript_disease_relevant:
                comment = (
                    f"Transcript {transcript.transcript_id} is not disease relevant."
                )
                result = False
            elif (
                transcript.diff_len_protein_percent <= threshold_diff_len_prot_percent
                and transcript.len_change_in_repetitive_region
            ):
                comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region overlaps repetitive region."
                result = True
            else:
                comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}."
                result = False
            rule_result = RuleResult(
                "BP3",
                rule_type.GENERAL,
                evidence_type.BENIGN,
                result,
                evidence_strength.SUPPORTING,
                comment,
            )
            results.append(rule_result)
        final_result = summarise_results_per_transcript(results)
        return final_result
