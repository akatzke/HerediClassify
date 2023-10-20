#!/usr/bin/env python3

from enum import Enum
import warnings

from refactoring.variant import PredictionTools


class Prediction_result(Enum):
    UNKNOWN = 0.5
    PATHOGENIC = 1
    BENIGN = 0


def assess_prediction_tool(prediction: PredictionTools, threshold: list[float]):
    if len(threshold) == 2:
        return assess_two_thresholds(prediction.value, threshold)
    elif len(threshold) == 1:
        return assess_one_threshold(prediction.value, threshold[0])
    else:
        warnings.warn(f"Number of thresholds supplied for {PredictionTools} wrong")


def assess_two_thresholds(data: float, thresholds: list[float]) -> Prediction_result:
    """ """
    if data <= thresholds[0]:
        return Prediction_result.BENIGN
    elif data >= thresholds[1]:
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
