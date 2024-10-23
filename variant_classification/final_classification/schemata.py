#!/usr/bin/env python3

import rule_combinations as rc

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

schema_atm = schema_acmg
schema_brca1 = schema_acmg
schema_brca2 = schema_acmg
schema_cdh1 = schema_acmg
schema_palb2 = schema_acmg
schema_pten = schema_acmg
schema_tp53 = schema_acmg
