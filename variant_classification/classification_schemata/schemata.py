#!/usr/bin/env python3

import classification_schemata.rule_combinations as rc

schema_acmg = {
    1: [rc.benign_1, rc.benign_2],
    2: [rc.likely_benign_1, rc.likely_benign_2],
    4: [
        rc.likely_pathogenic_1,
        rc.likely_pathogenic_2,
        rc.likely_pathogenic_3,
        rc.likely_pathogenic_4,
        rc.likely_pathogenic_5,
        rc.likely_pathogenic_6,
        rc.likely_pathogenic_7,
    ],
    5: [
        rc.pathogenic_1,
        rc.pathogenic_2,
        rc.pathogenic_3,
        rc.pathogenic_4,
        rc.pathogenic_5,
        rc.pathogenic_6,
        rc.pathogenic_7,
        rc.pathogenic_8,
        rc.pathogenic_9,
    ],
}

schema_atm = schema_acmg.copy()
schema_atm.update({2: [rc.likely_benign_1, rc.likely_benign_2, rc.likely_benign_atm_1]})

schema_brca1 = {
    1: [
        rc.benign_1,
        rc.benign_2,
        rc.benign_brca1_1,
        rc.benign_brca1_2,
        rc.benign_brca1_3,
    ],
    2: [
        rc.likely_benign_1,
        rc.likely_benign_2,
        rc.likely_benign_brca1_1,
        rc.likely_benign_brca1_2,
        rc.likely_benign_brca1_3,
    ],
    4: [
        rc.pathogenic_6,
        rc.likely_pathogenic_brca1_1,
        rc.likely_pathogenic_brca1_2,
        rc.likely_pathogenic_3,
        rc.likely_pathogenic_6,
        rc.likely_pathogenic_7,
    ],
    5: [
        rc.pathogenic_1,
        rc.pathogenic_2,
        rc.pathogenic_brca1_1,
        rc.pathogenic_5,
        rc.pathogenic_brca1_2,
        rc.pathogenic_brca1_3,
        rc.pathogenic_brca1_4,
        rc.pathogenic_7,
        rc.pathogenic_8,
        rc.pathogenic_9,
    ],
}

schema_brca2 = schema_brca1.copy()
schema_brca2.update(
    {
        4: [
            rc.likely_pathogenic_4,
            rc.likely_pathogenic_5,
            rc.likely_pathogenic_6,
            rc.likely_pathogenic_7,
        ]
    }
)


schema_cdh1 = schema_acmg.copy()

schema_palb2 = schema_atm.copy()

schema_pten = schema_atm.copy()

schema_tp53 = schema_acmg.copy()
