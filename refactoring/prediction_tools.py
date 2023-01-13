#!/usr/bin/env python3
from numpy import mean
from enum import Enum
import pandas as pd

from import_config import Prediction_tool_threshold


class Prediction_result(Enum):
    UNKNOWN = 0.5
    PATHOGENIC = 1
    BENIGN = 0


def get_splicing_prediction(
    data: pd.Series, threshold: Prediction_tool_threshold
) -> Prediction_result:
    """
    Get results from all given splicing predicition tools and summarize their results
    """
    spliceai_result = parse_spliceai_pred(data["SpliceAI"], threshold.SpliceAI)
    maxentscan_result = parse_maxentscan_pred(data["MaxEntScan"], threshold.MaxEntScan)
    return aggregate_splicing_predictions(spliceai_result, maxentscan_result)


def get_pathogenicity_prediction(
    data: pd.Series, threshold: Prediction_tool_threshold
) -> bool:
    """
    Get results from all given pathogenicity prediction tools and summarize their results
    """
    revel_result = assess_two_thresholds(
        data["REVEL"], threshold.revel_benign, threshold.revel_pathogenic
    )
    cadd_result = assess_one_threshold(data["CADD"], threshold.CADD)
    pyhlop_result = assess_one_threshold(data["phylop"], threshold.pyhlop)
    result = aggregate_predictions([revel_result, cadd_result, pyhlop_result])
    return result


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
        warn()
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
