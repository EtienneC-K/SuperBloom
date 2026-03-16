#!/usr/bin/env python3
from common import run_parameter_sweep_cli


if __name__ == "__main__":
    run_parameter_sweep_cli(
        parameter="n_hashes",
        default_values=[2, 4, 6, 8, 10, 12],
        description="Sweep number of hashes while keeping other parameters at defaults.",
    )
