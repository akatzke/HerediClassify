#!/usr/bin/env python3
import pandas as pd
from dataclasses import dataclass
from typing import Literal, Optional, Callable

from pandas.core.tools.datetimes import _return_parsed_timezone_results

from prediction_tools import get_prediction, Prediction_result

STRENGTH_TYPE = Literal["supporting", "moderate", "strong", "very strong"]


@dataclass
class Rule_result:
    name: str
    strength: STRENGTH_TYPE
    function: Callable
    comment: Optional[str] = None
    status: Optional[bool] = None

    def run(self, data: pd.Series):
        self.strength, self.comment, self.status = self.function(data)


def assess_pp3(
    self,
    data: pd.Series,
    threshold: pd.Series,
    splicing_tools: list,
    pathogenicity_tools: list,
) -> None:
    """
    Function assessing pp3: compuational evidence for variant being pathogenic
    """
    splicing_prediction = get_prediction(data, threshold, splicing_tools)
    pathogenicity_prediction = get_prediction(data, threshold, pathogenicity_tools)
    self.status = apply_logic_pp3(splicing_prediction, pathogenicity_prediction)
    self.commet = create_comment_prediction_tool(
        self.status, splicing_prediction, pathogenicity_prediction
    )


def create_comment_prediction_tool(
    self,
    splicing_prediction: Prediction_result,
    pathogenicity_prediction: Prediction_result,
) -> str:
    comment = f"Rule {self.name} has been assessed as {self.status}. "
    if splicing_prediction == Prediction_result.PATHOGENIC:
        comment = comment + "Splcing prediction predicts effect on splicing. "
    else:
        comment = comment + "Splcing prediction predicts no effect on splicing. "
    if pathogenicity_prediction == Prediction_result.PATHOGENIC:
        comment = comment + "Pathogenicity prediction predicts effect. "
    else:
        comment = comment + "Pathogenicity prediction predicts no effect. "
    return comment


def apply_logic_pp3(
    splicing_prediction: Prediction_result, pathogenicity_prediction: Prediction_result
) -> bool:
    if (
        splicing_prediction == Prediction_result.PATHOGENIC
        and pathogenicity_prediction == Prediction_result.PATHOGENIC
    ):
        return True
    elif (
        splicing_prediction == Prediction_result.BENIGN
        and pathogenicity_prediction == Prediction_result.BENIGN
    ):
        return False
    else:
        return False


def assess_bp4(
    self,
    data: pd.Series,
    thresholds: pd.Series,
    splicing_tools: list,
    pathogenicity_tools: list,
) -> None:
    """
    Function assessing BP4: computational evidence for variant being bening
    """
    splicing_prediction = get_prediction(data, thresholds, splicing_tools)
    pathogenicity_prediction = get_prediction(data, thresholds, pathogenicity_tools)
    self.status = apply_logic_bp4(splicing_prediction, pathogenicity_prediction)
    self.comment = create_comment_prediction_tool(
        self, splicing_prediction, pathogenicity_prediction
    )


def apply_logic_bp4(
    splicing_prediction: Prediction_result, pathogenicity_prediction: Prediction_result
) -> bool:
    if (
        splicing_prediction == Prediction_result.PATHOGENIC
        and pathogenicity_prediction == Prediction_result.PATHOGENIC
    ):
        return True
    elif (
        splicing_prediction == Prediction_result.BENIGN
        and pathogenicity_prediction == Prediction_result.BENIGN
    ):
        return False
    else:
        return False


def assess_bp7(
    self, data: pd.Series, thresholds: pd.Series, list_of_tools: list
) -> None:
    """
    Function assessing BP7: computational evidence for splicing in synonymous, missense variants
    """
    splicing_prediction = get_prediction(data, thresholds, list_of_tools)
    self.status = apply_logic_bp7(splicing_prediction)
    self.comment = f"Rule {self.name} "


def apply_logic_bp7(splicing_prediction: Prediction_result) -> bool:
    if splicing_prediction == Prediction_result.BENIGN:
        return True
    else:
        return False
