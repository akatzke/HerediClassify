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
    ## Setup status PVS1 and BA1
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
    ## Disable PP3 in case PVS1 applies
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
    ## Incompatibility BA1 and BS1
    if ba1_applies and rules.get("BS1", {}).get("status", False):
        rules["BS1"]["status"] = False
        new_comment = rules["BS1"]["comment"] + " BS1 does not apply when BA1 applies."
        rules["BS1"]["comment"] = new_comment

    ## ---- Start gene specific incompatibilities ----

    ## For ATM and PALB2
    if name_config == "ACMG ATM" or name_config == "ACMG PALB2":
        ## Incompatibility of PVS1_RNA and PM5_splicing
        if not rules.get("PVS1_RNA", {}).get("status", False) and rules.get(
            "PM5_splicing", {}
        ).get("status", False):
            rules["PM5_splicing"]["status"] = False
            new_comment = (
                rules["PM5_splicing"]["comment"]
                + " As PVS1_RNA does not apply to this variant, PM5 splicing can not be applied."
            )
            rules["PM5_splicing"]["comment"] = new_comment
        ## Incompatibility of BP7_RNA and BP4_splicing
        if rules.get("BP7_RNA", {}).get("status", False) and rules.get(
            "BP4_splicing", {}
        ).get("status", False):
            rules["BP4_splicing"]["status"] = False
            new_comment = (
                rules["BP4_splicing"]["comment"]
                + " As BP4 splicing can only apply if BP7_RNA does not apply, BP4 is set to False."
            )

    ## For BRCA1 and BRCA2
    if name_config == "ACMG BRCA1" or name_config == "ACMG BRCA2":
        ## Incomatibility of PVS1 and PM5
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
        ## Application of BP7 requries application of BP4
        if not rules.get("BP4_splicing", False) and rules.get("BP7_splicing", False):
            rules["BP7_splicing"]["status"] = False
            rules["BP7_splicing"]["comment"] = (
                rules["BP7_splicing"]["comment"]
                + " As BP7 can only apply if BP4 applies as well and BP4 does not apply, BP7 is set to False."
            )
    return rules
