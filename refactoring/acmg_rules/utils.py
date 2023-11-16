#!/usr/bin/env python3

from abc import ABC, abstractmethod
from typing import Callable
from dataclasses import dataclass
from enum import Enum

from refactoring.information import Classification_Info, Info


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


def summarise_results_per_transcript(results: list[RuleResult]) -> RuleResult:
    return results[0]
