"""Run the strict eight-case DSS gate on the biased-loop research model."""

from __future__ import annotations

import argparse
from dataclasses import asdict, fields
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

import RIVASandBaselineAudit as baseline  # noqa: E402
from RIVASandBiasedLoop import (  # noqa: E402
    RIVASandBiasedLoopModel,
    RIVASandBiasedLoopParameters,
)


strict = baseline.audit_module
strict.RIVASandDBFBModel = RIVASandBiasedLoopModel
strict.RIVASandDBFBParameters = RIVASandBiasedLoopParameters


def parameters_from_args(
    args: argparse.Namespace,
) -> RIVASandBiasedLoopParameters:
    defaults = asdict(RIVASandBiasedLoopParameters())
    allowed = {item.name for item in fields(RIVASandBiasedLoopParameters)}
    values = defaults.copy()
    values.update({
        "branch_compliance_gain": args.branch_gain,
        "branch_compliance_exponent": args.branch_exponent,
        "branch_directional_balance": args.branch_balance,
    })
    overrides = strict._parameter_overrides(args.parameter, defaults)
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise ValueError(f"unknown biased-loop parameter(s): {', '.join(unknown)}")
    values.update(overrides)
    return RIVASandBiasedLoopParameters(**values)


strict.parameters_from_args = parameters_from_args


def build_parser() -> argparse.ArgumentParser:
    parser = strict.build_parser()
    parser.description = __doc__
    parser.add_argument("--branch-gain", type=float, default=4.5)
    parser.add_argument("--branch-exponent", type=float, default=1.0)
    parser.add_argument("--branch-balance", type=float, default=0.03)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    result = strict.audit(args)
    print(strict.json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
