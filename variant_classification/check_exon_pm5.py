#!/usr/bin/env python3

import pathlib
import logging
from collections.abc import Callable
from typing import Optional

import pandas as pd

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)

from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo
from var_type import VARTYPE_GROUPS


logger = logging.getLogger("GenOtoScope_Classify.Check_exon_pm5")


def annotate_exon_classification_pm5(
    variant: VariantInfo,
    annotated_transcripts: list[TranscriptInfo],
    exon_pm5_path: pathlib.Path,
    mane_path: pathlib.Path,
) -> Optional[RuleResult]:
    """
    Check if variant is listed in the preclassified splice sites by VCEP
    Here the PM5 classification is accessed
    """
    results = {}
    for transcript in annotated_transcripts:
        if any(
            var_type in VARTYPE_GROUPS.EXONIC.value for var_type in transcript.var_type
        ):
            exon_pm5 = pd.read_csv(exon_pm5_path, sep="\t")
            exon_table_entries = exon_pm5[
                (exon_pm5.start <= transcript.ptc) & (exon_pm5.end >= transcript.ptc)
            ]
            if exon_table_entries.empty:
                return None
            elif exon_table_entries.shape[0] != 1:
                strongest_evidence_entry = select_entry_with_strongest_evidence(
                    exon_table_entries
                )
                result = RuleResult(
                    "PM5",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    bool(strongest_evidence_entry.rule_status),
                    evidence_strength(strongest_evidence_entry.evidence_strength),
                    str(strongest_evidence_entry.final_comment),
                )
            else:
                strongest_evidence_entry = exon_table_entries
                result = RuleResult(
                    "PM5",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    bool(strongest_evidence_entry.rule_status.values[0]),
                    evidence_strength(
                        strongest_evidence_entry.evidence_strength.values[0]
                    ),
                    strongest_evidence_entry.comment.values[0],
                )
            results[transcript.transcript_id] = result
    if len(results) == 0:
        if not annotated_transcripts:
            comment = (
                "No annotated transcripts provided, PM5_enigma can not be applied."
            )
        else:
            comment = f"Accessing PM5_enigma does not apply to this variant, as PM5_enigma does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
        result = RuleResult(
            "PM5",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.MODERATE,
            comment,
        )
    else:
        result = summarise_results_per_transcript(results, mane_path)
    return result


def get_annotate_exon_classification_pm5(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for checking for entries in splice site classification table
    """
    return (
        annotate_exon_classification_pm5,
        (
            class_info.VARIANT,
            class_info.ANNOTATED_TRANSCRIPT_LIST,
            class_info.EXON_PM5_PATH,
            class_info.MANE_TRANSCRIPT_LIST_PATH,
        ),
    )


def select_entry_with_strongest_evidence(data: pd.DataFrame) -> pd.Series:
    """
    From table select entry with strongest evidence strength
    In case all entries have the same evidence strength, return the first
    """
    if not data[data.evidence_strength == "very_strong"].empty:
        very_strong = data[data.evidence_strength == "very_strong"]
        if very_strong.shape[0] == 1:
            very_strong["final_comment"] = very_strong["comment"]
            return very_strong.iloc[0, :]
        else:
            very_strong_comment = summarise_comments(very_strong)
            return very_strong_comment.iloc[0, :]
    elif not data[data.evidence_strength == "strong"].empty:
        strong = data[data.evidence_strength == "strong"]
        if strong.shape[0] == 1:
            strong["final_comment"] = strong["comment"]
            return strong.iloc[0, :]
        else:
            strong_comment = summarise_comments(strong)
            return strong_comment.iloc[0, :]
    elif not data[data.evidence_strength == "moderate"].empty:
        moderate = data[data.evidence_strength == "moderate"]
        if moderate.shape[0] == 1:
            moderate["final_comment"] = moderate["comment"]
            return moderate.iloc[0, :]
        else:
            moderate_comment = summarise_comments(moderate)
            return moderate_comment.iloc[0, :]
    elif data[data.evidence_strength == "supporting"].empty:
        supporting = data[data.evidence_strength == "supporting"]
        if supporting.shape[0] == 1:
            supporting["final_comment"] = supporting["comment"]
            return supporting.iloc[0, :]
        else:
            supporting_comment = summarise_comments(supporting)
            return supporting_comment.iloc[0, :]
    else:
        raise ValueError(f"The evidence strength does not match.")


def summarise_comments(data: pd.DataFrame) -> pd.DataFrame:
    """
    In case of two entries with the same evidence strength
    Join both entries
    """
    final_comment = ""
    for _, entry in data.iterrows():
        final_comment = final_comment + " " + entry.comment
    data["final_comment"] = final_comment
    return data
