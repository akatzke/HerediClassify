#!/usr/bin/env python3
from typing import Tuple
from math import isnan
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
    spliceai_result = parse_spliceai_pred(
        data["SpliceAI"], Prediction_tool_threshold.SpliceAI
    )
    maxentscan_result = parse_maxentscan_pred()
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
    result = aggregate_patho_predictions([revel_result, cadd_result, pyhlop_result])
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


def aggregate_patho_predictions(patho_predictions: list) -> Prediction_result:
    """
    Aggregate pathogenicity score from different predictors
    to assess total pathogenicity
    Logic: average of all non negative pathogenicity scores,
    if positive sum => Pathogenic, else => Benign
    """
    no_negative_predictions = []
    for prediction in patho_predictions:
        if prediction > -1:
            no_negative_predictions.append(prediction)
    if len(no_negative_predictions) == 0:
        return Prediction_result.UNKNOWN
    elif mean(no_negative_predictions) >= 1:
        return Prediction_result.PATHOGENIC
    else:
        return Prediction_result.BENIGN


def parse_maxentscan_pred(data: pd.Series) -> int:
    """
    Parse MaxEntScan column for current variant row
    """
    if data["MaxEntScan"] and ">" in str(data["MaxEntScan"]):
        # if MaxEntScan not null and it is parse-able
        # calculate MaxEntScan ratio
        impacting_splicing_ratios, no_impacting_splicing_ratios = [], []
        for ref_obs_scores in str(data["MaxEntScan"]).split(","):
            ref_score, obs_score = ref_obs_scores.split(">")
            if ref_score.isdigit() and obs_score.isdigit():
                ref_obs_ratio = abs(
                    float(obs_score) - float(ref_score) / float(ref_score)
                )
            else:
                ref_obs_ratio = -1.0
            if ref_obs_ratio >= 0.15:
                impacting_splicing_ratios.append(str(ref_obs_ratio))
            else:
                no_impacting_splicing_ratios.append(str(ref_obs_ratio))
    else:
        # no MaxEntScan score found
        impacting_splicing_ratios = None
    return impacting_splicing_ratios


def parse_dbscsnv_pred(data: pd.Series, threshold=float) -> int:
    """
    Parse dbscSNV splicing predictor column of current variant row
    """
    impacting_splicing_scores = []
    if data["dbscSNV"] and "/" in str(data["dbscSNV"]):
        ada_score, rf_score = str(data["dbscSNV"]).split("/")
        if ada_score.isdigit() and float(ada_score) > 0.6:
            impacting_splicing_scores.append(float(ada_score))
        if rf_score.isdigit() and float(rf_score) > 0.6:
            impacting_splicing_scores.append(float(rf_score))
    else:
        # no dbscSNV score found
        impacting_splicing_scores = None
    return impacting_splicing_scores


def aggregate_splicing_predictions(splicing_predictions) -> Prediction_result:
    """
    Aggregate splicing predictions to assess splicing impact
    Logic: all predictors, with score, should have a impacting score
    """
    sum_impact = 0
    for predictor, impacting_splicing_scores in splicing_predictions.items():
        if impacting_splicing_scores is not None and len(impacting_splicing_scores) > 0:
            sum_impact = sum_impact + 1
        else:
            sum_impact = sum_impact - 1
    if sum_impact > 0:
        return Prediction_result.PATHOGENIC
    elif sum_impact:
        return Prediction_result.UNKNOWN
    else:
        return Prediction_result.BENIGN


def aggregate_prediction_results(results: list[Prediction_result]) -> Prediction_result:
    if all(result == Prediction_result.BENIGN for result in results):
        return Prediction_result.BENIGN
    elif all(result == Prediction_result.PATHOGENIC for result in results):
        return Prediction_result.PATHOGENIC
    else:
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


def prepare_hbond(str) -> dict:
    split_list = str.split(
        "|",
    )
    predictions = {"test": split_list}
    return predictions
