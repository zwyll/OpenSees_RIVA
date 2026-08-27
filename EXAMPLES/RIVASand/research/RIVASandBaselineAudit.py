"""Private eight-case DSS audit for the frozen RIVA-Sand baseline.

This adapter deliberately reuses only the strict case matrix, stress-control
driver, metrics, and acceptance gates developed during the rejected DBFB
study.  It replaces that study's constitutive class and parameter factory
with the canonical standalone RIVA-Sand oracle.  No rejected constitutive
mechanism is imported into the response.

The file is a research harness and is not public OpenSees documentation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, fields
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
STRICT_PATH = (
    WORKSPACE / "OpenSees-rivasand-geostatic" / "EXAMPLES" / "RIVASand"
    / "research" / "RIVASandDBFBAudit.py"
)
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(STRICT_PATH.parent))

from rivasand_port.model import RIVASandModel, RIVASandParameters  # noqa: E402
from rivasand_port.reference import reference_parameters  # noqa: E402


SPEC = importlib.util.spec_from_file_location(
    "_rivasand_strict_baseline_metrics", STRICT_PATH
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(f"cannot load strict DSS audit from {STRICT_PATH}")
audit_module = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = audit_module
SPEC.loader.exec_module(audit_module)

# These are the only model hooks used by the strict audit.
audit_module.RIVASandDBFBModel = RIVASandModel
audit_module.RIVASandDBFBParameters = RIVASandParameters


def parameters_from_args(args: argparse.Namespace) -> RIVASandParameters:
    """Build the frozen baseline plus explicit command-line overrides."""
    defaults = asdict(reference_parameters())
    allowed = {item.name for item in fields(RIVASandParameters)}
    overrides = audit_module._parameter_overrides(args.parameter, defaults)
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise ValueError(f"unknown RIVA-Sand parameter(s): {', '.join(unknown)}")
    return reference_parameters(**overrides)


audit_module.parameters_from_args = parameters_from_args


def build_parser() -> argparse.ArgumentParser:
    parser = audit_module.build_parser()
    parser.description = __doc__
    for action in parser._actions:
        if action.dest == "parameter":
            action.help = "override a frozen RIVA-Sand parameter (repeatable)"
    return parser


def main() -> None:
    args = build_parser().parse_args()
    result = audit_module.audit(args)
    print(audit_module.json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
