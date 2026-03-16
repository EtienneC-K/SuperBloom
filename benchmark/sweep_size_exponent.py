#!/usr/bin/env python3
from common import run_parameter_sweep_cli


if __name__ == "__main__":
    run_parameter_sweep_cli(
        parameter="size_exponent",
        default_values=[31, 33, 35],
        description="Sweep total filter memory exponent while keeping other parameters at defaults.",
    )
