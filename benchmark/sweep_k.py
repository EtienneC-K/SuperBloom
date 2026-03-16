#!/usr/bin/env python3
from common import run_parameter_sweep_cli


if __name__ == "__main__":
    run_parameter_sweep_cli(
        parameter="k",
        default_values=[21, 25, 31, 41, 51],
        description=(
            "Sweep k while keeping other parameters at defaults. "
            "s is auto-adjusted to k-4 for each point."
        ),
        link_s_to_k=True,
    )
