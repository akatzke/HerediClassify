#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)
from information import Info, Classification_Info
from variant import TranscriptInfo, VariantInfo
from transcript_annotated import TranscriptInfo_exonic, TranscriptInfo_intronic


class Pm5_protein_ptc(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    Implementing the rule specifications for ATM
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.ANNOTATED_TRANSCRIPT_LIST,
                class_info.POS_LAST_KNOWN_PATHO_PTC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        annotated_transcript: list[TranscriptInfo],
        pos_last_known_patho_ptc_dict: dict[str, int],
    ) -> RuleResult:
        results = []
        for transcript in annotated_transcript:
            if isinstance(transcript, TranscriptInfo_exonic):
                try:
                    pos_last_known_patho_ptc = pos_last_known_patho_ptc_dict[
                        transcript.transcript_id
                    ]
                except KeyError:
                    raise KeyError(
                        f"Transcript {transcript.transcript_id} not in disease relevant transcripts: {pos_last_known_patho_ptc_dict.keys()}. Transcript should have been filtered out earlier."
                    )
                if transcript.ptc <= pos_last_known_patho_ptc:
                    status = True
                    comment = f"PTC ({transcript.ptc}) caused by variant is located upstream of last known pathogenic PTC {pos_last_known_patho_ptc}."
                else:
                    status = False
                    comment = f"PTC ({transcript.ptc}) caused by variant is located downstream of last known pathogenic PTC {pos_last_known_patho_ptc}."
                result = RuleResult(
                    "PM5",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    status,
                    evidence_strength.SUPPORTING,
                    comment,
                )
                results.append(result)
        if len(results) == 0:
            comment = f"PM5 does not apply to this variant, as PM5 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            final_result = RuleResult(
                "PM5",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                comment,
            )
        else:
            final_result = summarise_results_per_transcript(results)
        return final_result


class Pm5_splicing_ptc(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    Implementing the rule specifications for ATM
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.ANNOTATED_TRANSCRIPT_LIST,
                class_info.POS_LAST_KNOWN_PATHO_PTC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        annotated_transcript: list[TranscriptInfo],
        pos_last_known_patho_ptc_dict: dict[str, int],
    ) -> RuleResult:
        results = []
        for transcript in annotated_transcript:
            if isinstance(transcript, TranscriptInfo_intronic):
                try:
                    pos_last_known_patho_ptc = pos_last_known_patho_ptc_dict[
                        transcript.transcript_id
                    ]
                except KeyError:
                    raise KeyError(
                        f"Transcript {transcript.transcript_id} not in disease relevant transcripts: {pos_last_known_patho_ptc_dict.keys()}. Transcript should have been filtered out earlier."
                    )
                if transcript.ptc <= pos_last_known_patho_ptc:
                    status = True
                    comment = f"PTC ({transcript.ptc}) caused by variant is located upstream of last known pathogenic PTC {pos_last_known_patho_ptc}."
                else:
                    status = False
                    comment = f"PTC ({transcript.ptc}) caused by variant is located downstream of last known pathogenic PTC {pos_last_known_patho_ptc}."
                result = RuleResult(
                    "PM5",
                    rule_type.SPLICING,
                    evidence_type.PATHOGENIC,
                    status,
                    evidence_strength.SUPPORTING,
                    comment,
                )
                results.append(result)
        if len(results) == 0:
            comment = f"PM5 does not apply to this variant, as PM5 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            final_result = RuleResult(
                "PM5",
                rule_type.SPLICING,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                comment,
            )
        else:
            final_result = summarise_results_per_transcript(results)
        return final_result
