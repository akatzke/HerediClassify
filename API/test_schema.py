#!/usr/bin/env python3

import pandas as pd

path_example_input = "/home/katzkean/variant_classification/API/example_input.json"
path_schema_input = "/home/katzkean/variant_classification/API/schema_input.json"

def load_json(path: str) -> pd.DataFrame:
    return pd.json_normalize(path)
