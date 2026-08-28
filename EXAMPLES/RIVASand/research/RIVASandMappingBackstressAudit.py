"""Run the strict eight-case DSS gate on the mapping-backstress successor."""

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
from RIVASandMappingBackstress import (  # noqa: E402
    RIVASandMappingBackstressModel,
    RIVASandMappingBackstressParameters,
)


strict = baseline.audit_module
strict.RIVASandDBFBModel = RIVASandMappingBackstressModel
strict.RIVASandDBFBParameters = RIVASandMappingBackstressParameters


def parameters_from_args(
    args: argparse.Namespace,
) -> RIVASandMappingBackstressParameters:
    defaults = asdict(RIVASandMappingBackstressParameters())
    allowed = {item.name for item in fields(RIVASandMappingBackstressParameters)}
    overrides = strict._parameter_overrides(args.parameter, defaults)
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise ValueError(
            f"unknown mapping-backstress parameter(s): {', '.join(unknown)}"
        )
    defaults.update(overrides)
    return RIVASandMappingBackstressParameters(**defaults)


strict.parameters_from_args = parameters_from_args


def build_parser() -> argparse.ArgumentParser:
    parser = strict.build_parser()
    parser.description = __doc__
    return parser


def main() -> None:
    args = build_parser().parse_args()
    result = strict.audit(args)
    print(strict.json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
