#!/usr/bin/env python3
from enum import Enum
import pandas as pd
import warnings

from import_config import Prediction_tool_threshold


class Prediction_result(Enum):
    UNKNOWN = 0.5
    PATHOGENIC = 1
    BENIGN = 0


def get_prediction(data: pd.Series, threshold: pd.Series, list_of_tools: list):
    results = []
    for tool in list_of_tools:
        if len(list(filter(lambda x: tool in x, list(threshold.index)))) == 1:
            results.append(
                assess_one_threshold(data=data[tool], threshold=threshold[tool])
            )
        elif len(list(filter(lambda x: tool in x, list(threshold.index)))) == 2:
            threshold_benign = get_threshold(threshold, tool, "benign")
            threshold_pathogenic = get_threshold(threshold, tool, "pathogenic")
            results.append(
                assess_two_thresholds(
                    data=data[tool],
                    threshold_benign=threshold_benign,
                    threshold_pathogenic=threshold_pathogenic,
                )
            )
        elif len(list(filter(lambda x: tool in x, list(threshold.index)))) == 0:
            warnings.warn(
                f"No threshold was found for {tool}. {tool} is skipped for analysis"
            )
        else:
            warnings.warn(
                f"Too many thresholds were found for {tool}. Check data import."
            )
    return aggregate_prediction_results(results)


def get_threshold(thresholds: pd.Series, tool: str, type: str) -> int:
    name_threshold = tool + "_" + type
    print(name_threshold)
    return thresholds[name_threshold]


def assess_two_thresholds(
    data: float, threshold_benign: float, threshold_pathogenic: float
) -> Prediction_result:
    """ """
    if data <= threshold_benign:
        return Prediction_result.BENIGN
    elif data >= threshold_pathogenic:
        return Prediction_result.PATHOGENIC
    else:
        return Prediction_result.UNKNOWN


def assess_one_threshold(data: float, threshold: float) -> Prediction_result:
    """
    Given one threshold, check if
    """
    if data >= threshold:
        return Prediction_result.PATHOGENIC
    elif data < threshold:
        return Prediction_result.BENIGN
    else:
        warnings.warn("Threshold could not be compared to data")
        return Prediction_result.UNKNOWN


def aggregate_prediction_results(results: list[Prediction_result]) -> Prediction_result:
    if all(result == Prediction_result.BENIGN for result in results):
        return Prediction_result.BENIGN
    elif all(result == Prediction_result.PATHOGENIC for result in results):
        return Prediction_result.PATHOGENIC
    else:
        return Prediction_result.UNKNOWN


def parse_dbscsnv_pred(dbscsnv: str, threshold: float) -> Prediction_result:
    """
    Parse dbscSNV splicing predictor column of current variant row
    """
    result = []
    if dbscsnv and "/" in str(dbscsnv):
        ada_score, rf_score = str(dbscsnv).split("/")
        if ada_score.isdigit() and float(ada_score) > threshold:
            result.append(Prediction_result.PATHOGENIC)
        if rf_score.isdigit() and float(rf_score) > threshold:
            result.append(Prediction_result.PATHOGENIC)
    else:
        result.append(Prediction_result.UNKNOWN)
    if not result:
        return Prediction_result.BENIGN
    elif len(result) == 1:
        return result[0]
    else:
        if sum(result) == 2:
            return Prediction_result.PATHOGENIC
        elif sum(result) == 1:
            return Prediction_result.UNKNOWN


def parse_spliceai_pred(spliceai: str, threshold: float) -> Prediction_result:
    result = []
    single_scores = spliceai.split("|")[2:6]
    for score in single_scores:
        result.append(assess_one_threshold(float(score), threshold.SpliceAI))
    if sum(result) == 0:
        return Prediction_result.BENIGN
    elif sum(result) == 1:
        return Prediction_result.PATHOGENIC
    elif sum(result) > 1:
        warn(
            "Two or more prediction from SpliceAI meet the threshold criterium. Please manually check."
        )
        return Prediction_result.UNKNOWN
    else:
        warn("Something went very wrong in SpliceAI evaluation")
        return Prediction_result.UNKNOWN
