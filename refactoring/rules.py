#!/usr/bin/env python3
import pandas as pd
from dataclasses import dataclass
from typing import Literal, Optional, Callable

from prediction_tools import get_pathogenicity_prediction, get_splicing_prediction

STRENGTH_TYPE = Literal["supporting", "moderate", "strong", "very strong"]


@dataclass
class Rule_result:
    name: str
    strength: STRENGTH_TYPE
    comment: Optional[str]
    status: Optional[bool]
    function: Callable

    def run(self, data: pd.Series):
        self.strength, self.comment, self.status = self.function(data)


class Assess_pp3(rule):
    name: str = "PP3"
    strength: STRENGTH_TYPE = "supporting"
    comment: Optional[str] = None
    status: Optional[bool] = None


def assess_pp3(self, data: pd.Series) -> None:
    """
    Function assessing pp3: compuational evidenve for variant being pathogenic
    """
    splicing_prediction = get_splicing_prediction(data)
    pathogenicity_prediction = get_pathogenicity_prediction(data)
    self.status = splicing_prediction or pathogenicity_prediction
    self.comment = ""


def assess_bp4(data: pd.Series) -> dict:
    """
    Function assessing BP4: computational evidence for variant being bening
    """
    splicing_prediction = get_splicing_prediction(data)
    pathogenicity_prediction = get_pathogenicity_prediction(data)
    does_bp4_apply = splicing_prediction and pathogenicity_prediction
    comment_bp4 = ""
    return {"rules_status": does_bp4_apply, "rules_comment": comment_bp4}


def assess_bp7(self, data: pd.Series) -> None:
    """
    Function assessing BP7: computational evidence for splicing in synonymous, missense variants
    """
    splicing_prediction = get_splicing_prediction(data)
    self.status = not splicing_prediction
    self.comment = ""
