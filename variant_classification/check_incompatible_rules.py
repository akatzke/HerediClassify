#!/usr/bin/env python3

from acmg_rules.utils import evidence_strength, rule_type


def check_incompatible_rules(
    rules: dict[str, dict], name_config: str, rules_list: list[str]
) -> dict[str, dict]:
    """
    Filter out incompatible rules
    E.g. {PVS1 is incompatible with PM4
    PVS1_splicing is incompatible with PS1_splicing
    """
    ba1_applies = rules.get("BA1", {}).get("status", False)
    pvs1_applies = rules.get("PVS1_splicing", {}).get("status", False) or rules.get(
        "PVS1_protein", {}
    ).get("status", False)
    pvs1_splicing = rules.get("PVS1_splicing", {}).get("status", False)
    pvs1_splicing_very_strong = (
        pvs1_splicing
        and rules.get("PVS1", {}).get("strength", None)
        == evidence_strength.VERY_STRONG.value
    )
    ## Disable PM4 in case PVS1 applies
    if pvs1_applies and rules.get("PM4", {}).get("status", False):
        rules["PM4"]["status"] = False
        new_comment = (
            rules["PM4"]["comment"]
            + " PM4 does not apply, as PVS1 already applies to the variant."
        )
        rules["PM4"]["comment"] = new_comment
    ## In case ClinGen specifications for splicing apply
    ## Adjust evidence strength of PS1 in case PVS1 splicing very strong applies
    if (
        pvs1_splicing_very_strong
        and "ps1_splicing_clingen" in rules_list
        and (
            rules.get("PS1_splicing", {}).get("status", False)
            and rules.get("PS1_splicing", {}).get("rule_type", None)
            == rule_type.SPLICING.value
        )
    ):
        rules["PS1_splicing"]["strength"] = evidence_strength.SUPPORTING.value
        new_comment = (
            rules["PS1_splicing"]["comment"]
            + " Correct evidence strength to supporting as PVS1 splicing applies with very strong evidence strength."
        )
        rules["PS1_splicing"]["coment"] = new_comment
    if name_config == "ACMG ATM":
        if pvs1_splicing and (
            rules.get("PP3_splicing", {}).get("status", False)
            and rules.get("PP3_splicing", {}).get("rule_type", None)
            == rule_type.SPLICING.value
        ):
            rules["PP3_splicing"]["status"] = False
            new_comment = (
                rules["PP3_splicing"]["comment"]
                + " PP3 splicing does not apply, as PVS1 splicing already applies to the variant."
            )
            rules["PP3_splicing"]["comment"] = new_comment
        if (
            rules.get("BP7_splicing", {}).get("status", False)
            and rules.get("BP7_splicing", {}).get("strength", None)
            != evidence_strength.SUPPORTING.value
            and rules.get("BP4_splicing", {}).get("status", False)
        ):
            rules["BP4_splicing"]["status"] = False
            new_comment = (
                rules["BP4_splicing"]["comment"]
                + " As BP7 RNA analysis applies to this variant. Do not apply BP4 splicing."
            )
            rules["BP4_splicing"]["comment"] = new_comment
    if name_config == "ACMG BRCA1" or name_config == "ACMG BRCA2":
        if (
            not pvs1_applies
            and rules.get("PVS1", {}).get("rule_type") != rule_type.SPLICING.value
            and rules.get("PM5", {}).get("status", False)
        ):
            rules["PM5"]["status"] = False
            new_comment = (
                rules["PM5"]["comment"]
                + " As PM5 can only apply if PVS1 applies as well and PVS1 does not apply, PM5 is set to False."
            )
            rules["PM5"]["comment"] = new_comment
    if ba1_applies and rules.get("BS1", {}).get("status", False):
        rules["BS1"]["status"] = False
        new_comment = rules["BS1"]["comment"] + " BS1 does not apply when BA1 applies."
        rules["BS1"]["comment"] = new_comment
    return rules
