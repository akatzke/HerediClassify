#!/usr/bin/env python3

import pathlib
import pandas as pd


def format_VCEP_splice_table(path: pathlib.Path) -> None:
    """
    Excpets a csv formatted as follows
    - position: "c.20+1G>"
    - alternative_allele: "A, C, T"
    - rule_status: "FALSE"
    - evidence_strength: strong
    - comment: "Some comment"
    And returns the following:
    - position: "20+1"
    - alternative_allele: "A"
    - rule_status: "FALSE"
    - evidence_strength: strong
    - comment: "Some comment"
    """
    splice_table = pd.read_csv(path)
    ## Reformat position from "c.20+1G>" to "20+1"
    splice_table.position = splice_table.position.str.replace("c.", "")
    splice_table.position = splice_table.position.str.replace(" ", "")
    splice_table.position = splice_table.position.str[:-2]
    ## Reformat alternative_allele
    splice_table.alternative_allele = splice_table.alternative_allele.str.replace(
        " ", ""
    )
    splice_table.alternative_allele = splice_table.alternative_allele.str.split(",")
    splice_table = splice_table.explode("alternative_allele")
    ## Create path to reformatted table
    file_name = path.stem + "_reformatted.csv"
    out_path = path.parent / file_name
    splice_table.to_csv(out_path, index=False, sep="\t")
