#!/usr/bin/env python3

import pathlib

import pandas as pd


def format_blosum62(path: pathlib.Path) -> None:
    """
    Format blosum62 score
    Downloaded from: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
    """
    blosum = pd.read_csv(path)
    blosum_melt = blosum.melt(var_name="alt_aa", value_name="score", ignore_index=False)
    blosum_melt["ref_aa"] = blosum_melt.index
    blosum_out = blosum_melt[["ref_aa", "alt_aa", "score"]]
    out_name = path.stem + "_formatted.csv"
    out_path = path.parent / out_name
    blosum_out.to_csv(out_path, index=False, sep="\t")


def format_grantham(path: pathlib.Path) -> None:
    """
    Format grantham score
    Downloaded from: https://gist.github.com/danielecook/501f03650bca6a3db31ff3af2d413d2a
    """
    grantham = pd.read_csv(path, sep="\t")
    grantham_index = grantham.set_index("FIRST")
    # Add missing row and column
    grantham_index["S"] = 0
    grantham_transposed = grantham_index.transpose()
    grantham_transposed["W"] = 0
    grantham_symmetrical = grantham_transposed + grantham_transposed.T
    #
    grantham_melt = grantham_index.melt(
        var_name="alt_aa", value_name="score", ignore_index=False
    )
    grantham_melt["ref_aa"] = grantham_melt.index
    grantham_out = grantham_melt[["ref_aa", "alt_aa", "score"]]
    out_name = path.stem + "_formatted.csv"
    out_path = path.parent / out_name
    grantham_out.to_csv(out_path, index=False, sep="\t")
