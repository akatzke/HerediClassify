#!/usr/bin/env python3
from rules import *
import pandas as pd


def test_fail_pp3() -> None:
    path_test_data = "/home/katzkean/mhh-lw-o/Variant_Classification/Test_files/prediction_test_fail.csv"
    test_data = pd.read_csv(path_test_data, index_col=False, header=0)
    test_data = pd.Series(
        test_data.iloc[
            0,
        ]
    )
    assert not get_pathogenicity_prediction(test_data)


def test_return_true_pp3() -> None:
    path_test_data = "/home/katzkean/mhh-lw-o/Variant_Classification/Test_files/prediction_test_true.csv"
    test_data = pd.read_csv(path_test_data, index_col=False, header=0)
    test_data = pd.Series(
        test_data.iloc[
            0,
        ]
    assert get_pathogenicity_prediction(test_data)


def test_return_false_pp3() -> None:
    path_test_data = "/home/katzkean/mhh-lw-o/Variant_Classification/Test_files/prediction_test_false.csv"
    test_data = pd.read_csv(path_test_data, index_col=False, header=0)
    test_data = pd.Series(
        test_data.iloc[
            0,
        ]
    assert not get_pathogenicity_prediction(test_data)
