#!/usr/bin/env python3
import pandas as pd
from utils import (
    parse_revel_pred,
    parse_cadd_pred,
    is_variant_conserved,
    parse_maxentscan_pred,
    parse_dbscsnv_pred,
)


def assess_PP3(data: pd.DataFrame) -> pd.DataFrame:
    """
    Function assessing PP3: computational evidence for pathogenicity
    """
    splicing_prediction = get_splicing_prediction()
    pathogenicity_prediction = get_pathogenicity_prediction()


def assess_BP4(data: pd.DataFrame) -> pd.DataFrame:
    """
    Function assessing BP4: computational evidence for variant being bening
    """
    splicing_prediction = get_splicing_prediction()
    pathogenicity_prediction = get_pathogenicity_prediction()


def assess_BP7(data: pd.DataFrame) -> pd.DataFrame:
    """
    Function assessing BP7: computational evidence for splicing in synonymous, missense variants
    """
    splicing_prediction = get_splicing_prediction()


def get_splicing_prediction(data: pd.DataFrame) -> pd.DataFrame:
    """
    Get results from all given splicing predicition tools and summarize their results
    """
    revel_result = parse_revel_pred()
    conserved_result = parse_pyhod


def get_pathogenicity_prediction(data: pd.DataFrame) -> pd.DataFrame:
    """
    Get results from all given pathogenicity prediction tools and summarize their results
    """
