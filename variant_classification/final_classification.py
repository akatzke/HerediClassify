#!/usr/bin/env python3

import pandas as pd

import classification_schemata.schemata as Class_schema
from classification_schemata.utils import (
    get_classifications_from_rule_combinations,
    get_final_classification_from_possible_classes,
)


def get_final_classifications(rules: dict, config: dict) -> dict:
    """
    Get final classification for variants
    """
    rules_df = pd.DataFrame(rules).transpose()
    applicable_rules = rules_df[rules_df.status == True]
    # Get final classification splice evidence
    rules_splicing = applicable_rules[
        applicable_rules.rule_type.isin(["splicing", "general"])
    ]
    class_splicing = get_classification(
        rules_splicing, config["name"], config["version"]
    )
    # Get final classification protein evidence
    rules_protein = applicable_rules[
        applicable_rules.rule_type.isin(["protein", "general"])
    ]
    class_protein = get_classification(rules_protein, config["name"], config["version"])
    # Add results to dictionary
    rules["classification protein"] = class_protein
    rules["classification_splicing"] = class_splicing
    return rules


def get_classification(rule_results: pd.DataFrame, config: str, version: str) -> int:
    """
    Execute final classification
    """
    schema = VERSION_CLASS_SCHEMATA.get(config, {}).get(version, None)
    if schema is None:
        raise ValueError(
            f"No final classification schemata defined for configuration {config} version {version}. Please check."
        )
    abs_rule_results = create_evidence_strength_count(rule_results)
    abs_rule_results["applicable_rules"] = rule_results.index
    possible_classes = get_classifications_from_rule_combinations(
        schema, abs_rule_results
    )
    final_class = get_final_classification_from_possible_classes(possible_classes)
    return final_class


def create_evidence_strength_count(rules: pd.DataFrame) -> dict:
    """
    From rules table create dictionary counting how often each evidence strenght is applicable
    """
    rules["evidence"] = rules.evidence_type + "_" + rules.strength
    counts = rules.groupby("evidence").count()
    return counts.rule_type.to_dict()


VERSION_CLASS_SCHEMATA = {
    "ACMG standard + SVI": {"1.0.0": Class_schema.schema_acmg},
    "ACMG ATM": {"1.3.0": Class_schema.schema_atm},
    "ACMG BRCA1": {"1.1.0": Class_schema.schema_brca1},
    "ACMG BRCA2": {"1.1.0": Class_schema.schema_brca2},
    "ACMG CDH1": {"3.1.0": Class_schema.schema_cdh1},
    "ACMG PALB2": {"1.1.0": Class_schema.schema_palb2},
    "ACMG PTEN": {"3.1.0": Class_schema.schema_pten},
    "ACMG TP53": {"1.4.0": Class_schema.schema_tp53},
}
