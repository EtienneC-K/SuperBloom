#!/usr/bin/env python3
from common import run_parameter_sweep_cli


if __name__ == "__main__":
    run_parameter_sweep_cli(
        parameter="threads",
        default_values=[1, 2, 4, 8, 16, 32],
        description="Sweep thread count while keeping algorithm parameters at defaults.",
    )
