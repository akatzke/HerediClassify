#!/usr/bin/env python3
import pandas as pd
from typing import Tuple
from math import isnan
from numpy import mean


def assess_PP3(data: pd.Series) -> dict:
    """
    Function assessing PP3: computational evidence for pathogenicity
    """
    splicing_prediction = get_splicing_prediction(data)
    pathogenicity_prediction = get_pathogenicity_prediction(data)
    does_pp3_apply = splicing_prediction or pathogenicity_prediction
    comment_pp3 = ""
    return {"rule_status": does_pp3_apply, "rule_comment": comment_pp3}


def assess_BP4(data: pd.Series) -> dict:
    """
    Function assessing BP4: computational evidence for variant being bening
    """
    splicing_prediction = get_splicing_prediction(data)
    pathogenicity_prediction = get_pathogenicity_prediction(data)
    does_bp4_apply = splicing_prediction and pathogenicity_prediction
    comment_bp4 = ""
    return {"rules_status": does_bp4_apply, "rules_comment": comment_bp4}


def assess_BP7(data: pd.Series) -> dict:
    """
    Function assessing BP7: computational evidence for splicing in synonymous, missense variants
    """
    splicing_prediction = get_splicing_prediction(data)
    does_bp7_apply = not splicing_prediction
    comment_bp7 = ""
    return {"rules_status": does_bp7_apply, "rules_comment": comment_bp7}


def get_splicing_prediction(data: pd.Series) -> Tuple[bool, str]:
    """
    Get results from all given splicing predicition tools and summarize their results
    """
    result = True
    comment = ""
    return result, comment


def get_pathogenicity_prediction(data: pd.Series) -> bool:
    """
    Get results from all given pathogenicity prediction tools and summarize their results
    """
    revel_result = parse_revel_pred(data)
    cadd_result = parse_cadd_pred(data)
    pyhlop_result = parse_pyhlop_pred(data)
    result = aggregate_patho_predictions([revel_result, cadd_result, pyhlop_result])
    return result


def parse_revel_pred(data: pd.Series) -> float:
    """
    Parse REVEL score from variant row
    and return its pathogenicity prediction
    """
    if isnan(data["REVEL"]):
        return -1
    elif float(data["REVEL"]) <= 0.15:
        return 0
    elif float(data["REVEL"]) >= 0.7:
        return 1
    else:
        return 0.5


def parse_pyhlop_pred(data: pd.Series) -> int:
    """
    Examine if variant is conserved using the PhyloP score
    and the cutoff of 1.6
    """
    if isnan(data["phyloP"]):
        return -1
    elif float(data["phyloP"]) > 1.6:
        return 1
    else:
        return 0


def parse_cadd_pred(data: pd.Series) -> int:
    """
    Parse CADD score from variant row
    and return its pathogenicity prediction
    """
    if isnan(data["CADD"]):
        return -1
    elif float(data["CADD"]) > 20.0:
        return 1
    else:
        return 0


def aggregate_patho_predictions(patho_predictions: list) -> str:
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
        return "unknown"
    elif mean(no_negative_predictions) >= 1:
        return "pathogenic"
    else:
        return "benign"


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


def parse_dbscsnv_pred(data: pd.Series) -> int:
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


def aggregate_splicing_predictions(splicing_predictions) -> bool:
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
        return True
    else:
        return False


def create_comment_computational_prediction(
    does_rule_apply: bool,
    pathogenicity_prediction: bool,
    splicing_predicition: bool,
    rule: str,
) -> str:

    return comment
