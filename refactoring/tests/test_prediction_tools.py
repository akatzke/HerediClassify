#!/usr/bin/env python3
from prediction_tools import (
    Prediction_result,
    assess_one_threshold,
    assess_two_thresholds,
    aggregate_prediction_results,
    get_prediction,
    get_threshold,
)
import pandas as pd

## Testing the basic function assessing whether or not the data meets the threshold


def test_one_threshold_pathogenic() -> None:
    assert assess_one_threshold(data=1.2, threshold=1.1) == Prediction_result.PATHOGENIC


def test_one_threshold_benign() -> None:
    assert assess_one_threshold(data=2.4, threshold=2.6) == Prediction_result.BENIGN


def test_one_threshold_equals_data_input() -> None:
    assert assess_one_threshold(data=1, threshold=1) == Prediction_result.PATHOGENIC


def test_two_thresholds_pathogenic() -> None:
    assert (
        assess_two_thresholds(data=2.2, threshold_benign=1.5, threshold_pathogenic=2.0)
        == Prediction_result.PATHOGENIC
    )


def test_two_thresholds_benign() -> None:
    assert (
        assess_two_thresholds(data=1.2, threshold_benign=1.5, threshold_pathogenic=2.0)
        == Prediction_result.BENIGN
    )


def test_two_thresholds_unknown() -> None:
    assert (
        assess_two_thresholds(data=1.7, threshold_benign=1.5, threshold_pathogenic=2.0)
        == Prediction_result.UNKNOWN
    )


## Testing the aggregation function


def test_aggregate_predictions_benign() -> None:
    results = [Prediction_result.BENIGN, Prediction_result.BENIGN]
    assert aggregate_prediction_results(results) == Prediction_result.BENIGN


def test_aggregate_predictions_pathogenic() -> None:
    results = [Prediction_result.PATHOGENIC, Prediction_result.PATHOGENIC]
    assert aggregate_prediction_results(results) == Prediction_result.PATHOGENIC


def test_aggregate_predictions_unknown() -> None:

    results = [
        Prediction_result.PATHOGENIC,
        Prediction_result.PATHOGENIC,
        Prediction_result.BENIGN,
    ]
    assert aggregate_prediction_results(results) == Prediction_result.UNKNOWN


## Testing the calling function


def test_get_predictions() -> None:
    thresholds = pd.Series(
        {"CADD": 24.0, "revel_benign": 1.2, "revel_pathogenic": 2.5, "phylop": 2.1}
    )
    tools = ["CADD", "revel", "phylop"]
    data = pd.Series({"CADD": 20.658, "revel": 2.6, "phylop": 1.8})
    assert (
        get_prediction(data=data, threshold=thresholds, list_of_tools=tools)
        == Prediction_result.UNKNOWN
    )
