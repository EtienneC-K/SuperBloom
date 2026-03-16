#!/usr/bin/env python3
from common import run_parameter_sweep_cli


if __name__ == "__main__":
    run_parameter_sweep_cli(
        parameter="s",
        default_values=[19, 23, 27, 29, 31],
        description="Sweep s-mer length s while keeping other parameters at defaults.",
    )
