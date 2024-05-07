#!/usr/bin/env python3

import pathlib

from abc import ABC, abstractmethod
from typing import Callable, Optional
from dataclasses import dataclass
from enum import Enum

import pandas as pd

from information import Classification_Info, Info


class evidence_strength(Enum):
    STAND_ALONE = "stand_alone"
    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"


class evidence_type(Enum):
    PATHOGENIC = "pathogenic"
    BENIGN = "benign"


class rule_type(Enum):
    GENERAL = "general"
    PROTEIN = "protein"
    SPLICING = "splicing"


@dataclass
class RuleResult:
    name: str
    type: rule_type
    evidence_type: evidence_type
    status: bool
    strength: evidence_strength
    comment: str

    def create_dict(self) -> dict:
        """
        Create a dictionary from the RuleResult object
        """
        rule_dict = {
            "rule_type": self.type.value,
            "evidence_type": self.evidence_type.value,
            "status": self.status,
            "strength": self.strength.value,
            "comment": self.comment,
        }
        rule_name = self.name.upper()
        if self.type is rule_type.PROTEIN or self.type is rule_type.SPLICING:
            rule_name = rule_name + f"_{self.type.value}"
        out_dict = {rule_name: rule_dict}
        return out_dict


class abstract_rule(ABC):
    @abstractmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        """
        Get function that assess rule
        """
        pass

    @abstractmethod
    def assess_rule(cls, args) -> RuleResult:
        """
        Assess rule
        """
        pass


def summarise_results_per_transcript(
    results: dict[str, RuleResult], mane_path: pathlib.Path
) -> RuleResult:
    # If dict has only one entry, return that entry as final result
    if len(results) == 1:
        return list(results.values())[0]
    # When there is mor
    mane_transcript_id = get_mane_transcript_id(list(results.keys()), mane_path)
    try:
        final_result = results[mane_transcript_id]
        del results[mane_transcript_id]
        compound_comment = construct_compound_comment_for_all_assessed_transcripts(
            results
        )
        final_result.comment = final_result.comment + " " + compound_comment
    except KeyError:
        compound_comment = construct_compound_comment_for_all_assessed_transcripts(
            results
        )
        final_result = RuleResult(
            "PVS1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.VERY_STRONG,
            comment=f"No MANE transcript has been classified for this variant. Therefore PVS1 is set to False. Other transcripts that were classified give the following assessments: "
            + compound_comment,
        )
    return final_result


def get_mane_transcript_id(
    transcript_ids: list[str], mane_path: pathlib.Path
) -> Optional[str]:
    """
    From a list of transcript, select get the MANE transcript id
    """
    mane_transcripts_df = pd.read_csv(mane_path, sep="\t")
    mane_transcripts = mane_transcripts_df.transcript.dropna()
    for transcript_id in transcript_ids:
        if any(mane_transcripts.str.contains(transcript_id)):
            return transcript_id
    return None


def construct_compound_comment_for_all_assessed_transcripts(
    results: dict[str, RuleResult]
) -> str:
    """
    From a dictionary of results, construct a compound comment
    """
    comment = ""
    for transcript_id, result in results.items():
        if result.status:
            tmp_comment = (
                transcript_id
                + ": "
                + result.name
                + " "
                + result.strength.value
                + " applies: "
                + result.comment
            )
        else:
            tmp_comment = (
                transcript_id
                + ": "
                + result.name
                + " does not applies: "
                + result.comment
            )
        comment = comment + tmp_comment + " "
    return comment
