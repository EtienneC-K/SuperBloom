#!/usr/bin/env python3
from common import run_parameter_sweep_cli


if __name__ == "__main__":
    run_parameter_sweep_cli(
        parameter="block_size_exponent",
        default_values=[8, 9, 10, 11, 12],
        description="Sweep block size exponent while keeping other parameters at defaults.",
    )
