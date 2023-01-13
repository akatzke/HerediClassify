#!/usr/bin/env python3
from prediction_tools import (
    assess_one_threshold,
    assess_two_thresholds,
    Prediction_result,
    aggregate_prediction_results,
)

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
