#!/usr/bin/env python3

from refactoring.rules import Assess_bp4, Assess_bp7, Assess_pp3, Rule_result


function_name_to_function = {"pp3": Assess_pp3, "bp4": Assess_bp4, "bp7": Assess_bp7}

if name not in function_name_to_function:
    raise ValueError("This rule could not be found!")

# Should be list of rule objects instead of list of dictionaries
my_rules = [
    {"name": "pp3", "strength": "low", "function": function_name_to_function["pp3"]}
]

for entry in my_rules:
    globals()[entry["name"] + "_result"] = Rule_result(name=entry[name])


def process_rules(
    names: list[str], strengths: list[str], rule_functions: list[Callable]
):
    for name, strength, rule_function in zip(names, strengths, rule_functions):
        for data_point in data:
            rule_function(data_point)
