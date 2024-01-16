#!/usr/bin/env python3

import pathlib
from collections.abc import Callable
from typing import Optional

import pandas as pd
from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
)

from information import Classification_Info, Info
from variant import TranscriptInfo
from var_type import VARTYPE_GROUPS


def annotate_exon_classification_pm5(
    transcripts: list[TranscriptInfo], path_exon_pm5: pathlib.Path
) -> Optional[RuleResult]:
    """
    Check if variant is listed in the preclassified splice sites by VCEP
    Here the PM5 classification is accessed
    """
    if len(transcripts) != 1:
        raise ValueError(
            "There should be only one disease relevant transcript defined for variants with a Splice Site PM5 Classificaton Table."
        )
    transcript = transcripts[0]
    if not any(
        var_type in VARTYPE_GROUPS.EXONIC.value for var_type in transcript.var_type
    ):
        return RuleResult(
            "PM5",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.VERY_STRONG,
            f"Accessing PM5_enigma does not apply to this variant, as PM5_enigma does not apply to variant types {', '.join([var_type.value for var_type in transcript.var_type])}.",
        )
    exon_pm5 = pd.read_csv(path_exon_pm5, sep="\t")
    exon_table_entries: pd.DataFrame = exon_pm5[
        (exon_pm5.start <= transcript.ptc) & (exon_pm5.end >= transcript.ptc)
    ]
    if exon_table_entries.empty:
        return None
    elif exon_table_entries.shape[0] != 1:
        strongest_evidence_entry = select_entry_with_strongest_evidence(
            exon_table_entries
        )
    else:
        strongest_evidence_entry = exon_table_entries
    return RuleResult(
        "PM5",
        rule_type.PROTEIN,
        evidence_type.PATHOGENIC,
        bool(strongest_evidence_entry.rule_status.values[0]),
        evidence_strength(strongest_evidence_entry.evidence_strength.values[0]),
        strongest_evidence_entry.comment.values[0],
    )


def get_annotate_exon_classification_pm5(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for checking for entries in splice site classification table
    """
    return (
        annotate_exon_classification_pm5,
        (class_info.ANNOTATED_TRANSCRIPT_LIST, class_info.EXON_PM5_PATH),
    )


def select_entry_with_strongest_evidence(data: pd.DataFrame) -> pd.DataFrame:
    """
    From table select entry with strongest evidence strength
    In case all entries have the same evidence strength, return the first
    """
    if not data[data.evidence_strength == "very_strong"].empty:
        very_strong = data[data.evidence_strength == "very_strong"]
        if very_strong.shape[0] == 1:
            return very_strong
        else:
            return very_strong.iloc[0, :]
    elif not data[data.evidence_strength == "strong"].empty:
        strong = data[data.evidence_strength == "strong"]
        if strong.shape[0] == 1:
            return strong
        else:
            return strong.iloc[0, :]
    elif not data[data.evidence_strength == "moderate"].empty:
        moderate = data[data.evidence_strength == "moderate"]
        if moderate.shape[0] == 1:
            return moderate
        else:
            return moderate.iloc[0, :]
    elif data[data.evidence_strength == "supporting"].empty:
        supporting = data[data.evidence_strength == "supporting"]
        if supporting.shape[0] == 1:
            return supporting
        else:
            return supporting.iloc[0, :]
    else:
        raise ValueError(f"The evidence strength does not match.")
