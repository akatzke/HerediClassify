#!/usr/bin/env python3

from acmg_rules.utils import RuleResult, evidence_strength, rule_type


def check_incompatible_rules(
    rules: list[RuleResult], name_config: str, rules_list: list[str]
) -> list[RuleResult]:
    """
    Filter out incompatible rules
    E.g. PVS1 is incompatible with PM4
    PVS1_splicing is incompatible with PS1_splicing
    """
    pvs1_applies, pvs1_splicing_very_strong_applies = False, False
    for rule in rules:
        if rule.name == "PVS1":
            pvs1_applies = True
            if (
                rule.evidence_type == rule_type.SPLICING
                and rule.strength == evidence_strength.VERY_STRONG
            ):
                pvs1_splicing_very_strong_applies = True
            break
    if pvs1_applies:
        for rule in rules:
            if rule.name == "PM4" and rule.status:
                rule.status = False
                rule.comment = (
                    rule.comment
                    + " The rule does not apply, as PVS1 already applies to the variant."
                )
    if pvs1_splicing_very_strong_applies and "ps1_splicing_clingen" in rules_list:
        for rule in rules:
            if rule.name == "PS1" and rule_type.SPLICING and rule.status:
                rule.strength = evidence_strength.SUPPORTING
                rule.comment = (
                    rule.comment
                    + " Correct evidence strength to supporting as PVS1 splicing applies with very strong evidence strength."
                )
    return rules
