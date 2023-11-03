#!/usr/bin/env python3

from abc import ABC, abstractmethod
from typing import Callable
from dataclasses import dataclass
from enum import Enum

import refactoring.information as info


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


class abstract_rule(ABC):
    @abstractmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
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


def summarise_results_per_transcript(results: list[RuleResult]) -> RuleResult:
    return results[0]
