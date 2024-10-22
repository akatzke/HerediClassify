#!/usr/bin/env python3

import pandas as pd

import final_classification as Class_schema


def get_final_classifications(rules: dict, config: dict) -> dict:
    """
    Get final classification for variants
    """
    rules_df = pd.DataFrame(rules).transpose()
    # Get final classification splice evidence
    rules_splicing = rules_df[rules_df.rule_type.isin(["splicing", "general"])]
    class_splicing = get_classification(
        rules_splicing, config["name"], config["version"]
    )
    # Get final classification protein evidence
    rules_protein = rules_df[rules_df.rule_type.isin(["protein", "general"])]
    class_protein = get_classification(rules_protein, config["name"], config["version"])


def get_classification(rules: pd.DataFrame, config: str, version: str) -> int:
    """
    Execute final classification
    """
    schemata = VERSION_CLASS_SCHEMATA.get(config, None).get(version, None)
    if schemata is None:
        raise ValueError(
            f"No final classification schemata defined for configuration {config} version {version}. Please check."
        )
    classification = schemata(rules)
    return classification


VERSION_CLASS_SCHEMATA = {
    "ACMG standard + SVI": {"1.0.0": Class_schema.final_classification_acmg},
    "ACMG ATM": {"1.3.0": Class_schema.final_classification_atm},
    "ACMG BRCA1": {"1.1.0": Class_schema.final_classification_brca1},
    "ACMG BRCA2": {"1.1.0": Class_schema.final_classification_brca2},
    "ACMG CDH1": {"3.1.0": Class_schema.final_classification_cdh1},
    "ACMG PALB2": {"1.1.0": Class_schema.final_classification_palb2},
    "ACMG PTEN": {"3.1.0": Class_schema.final_classification_pten},
    "ACMG TP53": {"1.4.0": Class_schema.final_classification_tp53},
}
