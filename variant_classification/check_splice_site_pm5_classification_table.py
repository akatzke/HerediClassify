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
)

from information import Classification_Info, Info
from variant import TranscriptInfo
from var_type import VARTYPE_GROUPS


def get_annotate_splice_site_classification_pm5(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for checking for entries in splice site classification table
    """
    return (
        annotate_splice_site_classification_pm5,
        (class_info.TRANSCRIPT, class_info.SPLICE_SITE_TABLE_PM5_PATH),
    )


logger = logging.getLogger("GenOtoScope_Classify.check_splice_table_pm5")


def annotate_splice_site_classification_pm5(
    transcripts: list[TranscriptInfo], path_splice_table: pathlib.Path
) -> Optional[RuleResult]:
    """
    Check if variant is listed in the preclassified splice sites by VCEP
    Here the PM5 classification is accessed
    """
    if len(transcripts) != 1:
        logger.warning(
            f"There should be only on disease relevant transcript defined for variants with a SpliceSite PM5 Classification Table."
        )
        return None
    transcript = transcripts[0]
    if not any(
        var_type in VARTYPE_GROUPS.INTRONIC.value for var_type in transcript.var_type
    ):
        return RuleResult(
            "PM5",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.VERY_STRONG,
            f"Accessing PM5_splice does not apply to this variant, as PM5_splice does not apply to variant types {', '.join([var_type.value for var_type in transcript.var_type])}.",
        )
    splice_table = pd.read_csv(path_splice_table, sep="\t")
    try:
        non_snp = ["ins", "del", "dup"]
        if any(entry in str(transcript.var_hgvs) for entry in non_snp):
            # In case of a larger deletion, find the position with the smallest offset from the exon
            if "_" in str(transcript.var_hgvs) and "del" in str(transcript.var_hgvs):
                if abs(transcript.var_hgvs.pos.start.offset) < abs(
                    transcript.var_hgvs.pos.end.offset
                ):
                    splice_table_entry = splice_table[
                        splice_table.position == str(transcript.var_hgvs.pos.start)
                    ]
                else:
                    splice_table_entry = splice_table[
                        splice_table.position == str(transcript.var_hgvs.pos.end)
                    ]

            # In case only a single nucleotide is deletet, inserted or duplicated
            else:
                splice_table_entry = splice_table[
                    splice_table.position == str(transcript.var_hgvs.pos)
                ]
        # Select correct row for single nucleotide exchanges
        else:
            splice_table_entry = splice_table[
                (splice_table.position == str(transcript.var_hgvs.pos))
                & (splice_table.alternative_allele == transcript.var_hgvs.edit.alt)
            ]
    except AttributeError:
        return None
    if splice_table_entry.empty:
        return None
    else:
        return RuleResult(
            "PM5",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            bool(splice_table_entry.rule_status.values[0]),
            evidence_strength(splice_table_entry.evidence_strength.values[0]),
            splice_table_entry.comment.values[0],
        )
