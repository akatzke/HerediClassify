#!/usr/bin/env python3

import re

import pandas as pd


def format_spliceai(df: pd.DataFrame) -> pd.DataFrame:
    df_nonan = df.dropna(subset=["SpliceAI"])
    df_nonan.SpliceAI = df_nonan.SpliceAI.str.replace("&", ",")
    two_spliceai = df_nonan[df_nonan.SpliceAI.str.contains(",", na=False)]
    one_spliceai = df_nonan[~df_nonan.SpliceAI.str.contains(",", na=False)]
    one_max_spliceai = format_spliceai_singel(one_spliceai)
    two_max_spliceai = format_spliceai_multiple(two_spliceai)
    df_spliceai_max = pd.concat([one_max_spliceai, two_max_spliceai])
    return df_spliceai_max


def format_spliceai_singel(df: pd.DataFrame) -> pd.DataFrame:
    """
    Format SpliceAI entry for a single splice site observation
    """
    df_non_empty = df[~df.SpliceAI.str.contains("\|\.\|", na=False)]
    splice_ai_entries = [
        "SpliceAI_observed",
        "SpliceAI_gene",
        "SpliceAI_score_1",
        "SpliceAI_score_2",
        "SpliceAI_score_3",
        "SpliceAI_score_4",
        "SpliceAI_pos_1",
        "SpliceAI_pos_2",
        "SpliceAI_pos_3",
        "SpliceAI_pos_4",
    ]
    df_non_empty[splice_ai_entries] = df_non_empty.SpliceAI.str.split("|", expand=True)
    df_non_empty["SpliceAI_max"] = df_non_empty[
        ["SpliceAI_score_1", "SpliceAI_score_2", "SpliceAI_score_3", "SpliceAI_score_4"]
    ].max(axis=1)
    df_non_empty.SpliceAI_max = pd.to_numeric(df_non_empty.SpliceAI_max)
    df_shortened = df_non_empty.drop(splice_ai_entries, axis=1)
    return df_shortened


def format_spliceai_multiple(df: pd.DataFrame) -> pd.DataFrame:
    """
    Format SpliceAI entry for multiple splice site observation
    """
    df["SpliceAI_max"] = 0
    for i, entry in df.iterrows():
        spliceai_list = entry.SpliceAI.split(",")
        local_max = []
        for spliceai in spliceai_list:
            if not re.match(spliceai, "\|\.\|"):
                scores = spliceai.split("|")[2:6]
                local_max.append(float(max(scores)))
        if local_max:
            global_max = max(local_max)
            df.loc[i, ["SpliceAI_max"]] = global_max
    return df
