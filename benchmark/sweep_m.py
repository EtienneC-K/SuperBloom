#!/usr/bin/env python3
from common import run_parameter_sweep_cli


if __name__ == "__main__":
    run_parameter_sweep_cli(
        parameter="m",
        default_values=[13, 17, 21, 25, 29],
        description="Sweep minimizer length m while keeping other parameters at defaults.",
    )
