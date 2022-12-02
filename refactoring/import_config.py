#!/usr/bin/env python3
import yaml

path_config = "/home/katzkean/variant_classification/config.yaml"


def load_config(path_yaml):
    with open(path_yaml) as config_file:
        config = yaml.safe_load(config_file)
    return config
