#!/usr/bin/env python3
from abc import ABC, abstractmethod
import pandas as pd
from dataclasses import dataclass
from typing import Literal, Optional, Callable

from prediction_tools import get_pathogenicity_prediction, get_splicing_prediction

STRENGTH_TYPE = Literal["supporting", "moderate", "strong", "very strong"]


class rule(ABC):
    @abstractmethod
    def __call__(self, data: pd.Series) -> dict:
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        pass

    @property
    @abstractmethod
    def strength(self) -> STRENGTH_TYPE:
        pass

    @property
    @abstractmethod
    def comment(self) -> str:
        pass

    @property
    @abstractmethod
    def status(self) -> bool:
        pass

@dataclass
class Rule_result:
    name:str
    strength: STRENGTH_TYPE
    comment: Optional[str]
    status: Optional[bool]
    function: Callable

    def run(self, data: pd.Series):
        self.strength, self.comment, self.status = self.function(data)

@dataclass
class rule2:
    name: str
    strength: STRENGTH_TYPE
    comment: Optional[str]
    status: Optional[bool]
    function: Callable[[pd.Series], [bool, str]]

    def run(self, data: pd.Series):
        self.status, self.comment = self.function(data)



@dataclass
class Assess_pp3(rule):
    name: str = "PP3"
    strength: STRENGTH_TYPE = "supporting"
    comment: Optional[str] = None
    status: Optional[bool] = None

    def __call__(self, data: pd.Series) -> None:
        """
        Function assessing PP3: computational evidence for pathogenicity
        """
        splicing_prediction = get_splicing_prediction(data)
        pathogenicity_prediction = get_pathogenicity_prediction(data)
        self.status = splicing_prediction or pathogenicity_prediction
        self.comment = ""


class Assess_bp4(rule):
    name: str = "BP4"
    strength: STRENGTH_TYPE = "supporting"
    comment: str = ""
    status: bool = False

    def __call__(data: pd.Series) -> dict:
        """
        Function assessing BP4: computational evidence for variant being bening
        """
        splicing_prediction = get_splicing_prediction(data)
        pathogenicity_prediction = get_pathogenicity_prediction(data)
        does_bp4_apply = splicing_prediction and pathogenicity_prediction
        comment_bp4 = ""
        return {"rules_status": does_bp4_apply, "rules_comment": comment_bp4}


class Assess_bp7(rule):
    name: str = "BP7"
    strength: STRENGTH_TYPE = "supporting"
    comment: str = ""
    status: bool = False

    def __call__(self, data: pd.Series) -> None:
        """
        Function assessing BP7: computational evidence for splicing in synonymous, missense variants
        """
        splicing_prediction = get_splicing_prediction(data)
        self.status = not splicing_prediction
        self.comment = ""


## Example
def example_implementation():
    data = pd.Series()
    assess_bp7 = Assess_bp7()
    assess_bp7(data)
    print (f"Rule {assess_bp7.name} returns {}")
