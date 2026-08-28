"""Untouched DSS holdout audit for the loose-biased shear-flow successor.

The original eight-case matrix was used repeatedly while the research
successor was developed.  This module therefore defines a separate, frozen
18-case matrix selected before looking at any successor predictions.  It
spans both DesignSafe projects, four density bands, zero through 0.375 static
bias, 40 and 100 kPa confinement, and CSR from 0.15 through 0.50.

The command can run the frozen production oracle, the signed-PT parent, or the
new successor through exactly the same stress-control driver and metrics.
This is private research infrastructure, not public OpenSees documentation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, fields
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

import RIVASandBaselineAudit as baseline  # noqa: E402
from RIVASandLooseBiasedShearFlow import (  # noqa: E402
    RIVASandLooseBiasedShearFlowModel,
    RIVASandLooseBiasedShearFlowParameters,
)
from RIVASandLooseUnbiasedCorrection import (  # noqa: E402
    RIVASandLooseUnbiasedCorrectionModel,
    RIVASandLooseUnbiasedCorrectionParameters,
)
from rivasand_port.model import RIVASandModel, RIVASandParameters  # noqa: E402
from rivasand_port.reference import reference_parameters  # noqa: E402


strict = baseline.audit_module
Case = strict.Case


# Experimental N values use the same 7.5%-strain definition as the original
# gate: double-amplitude strain for zero-bias tests and peak single-amplitude
# strain for biased tests.  Comparison cycles are fixed well before triggering.
HOLDOUT_CASES = (
    Case("h3484_dr55_u015", "PRJ-3484", "eo_0_631_sigv_40_CSR_0_150_Tau_0_.csv",
         40.0, 0.631, 0.15, 0.000, 17.0, "double_amplitude", 10.7188, 3),
    Case("h3484_dr71_u017", "PRJ-3484", "eo_0_588_sigv_40_CSR_0_170_Tau_0_.csv",
         40.0, 0.588, 0.17, 0.000, 26.0, "double_amplitude", 17.2266, 6),
    Case("h3484_dr76_u021", "PRJ-3484", "eo_0_576_sigv_40_CSR_0_210_Tau_0_.csv",
         40.0, 0.576, 0.21, 0.000, 16.0, "double_amplitude", 10.1563, 3),
    Case("h3484_p100_u017", "PRJ-3484", "eo_0_600_sigv_100_CSR_0_170_Tau_0_.csv",
         100.0, 0.600, 0.17, 0.000, 11.0, "double_amplitude", 7.15625, 2),
    Case("h3484_b025_c025", "PRJ-3484", "eo_0_601_sigv_40_CSR_0_250_Tau_10_.csv",
         40.0, 0.601, 0.25, 0.250, 22.0, "single_amplitude", 14.1797, 5),
    Case("h3484_b0375_c035", "PRJ-3484", "eo_0_601_sigv_40_CSR_0_350_Tau_15_.csv",
         40.0, 0.601, 0.35, 0.375, 22.0, "single_amplitude", 14.2031, 5),
    Case("h4666_dr50_u015", "PRJ-4666", "eo_0_658_sigv_40_CSR_0_150_Tau_0_.csv",
         40.0, 0.658, 0.15, 0.000, 12.0, "double_amplitude", 7.71094, 2),
    Case("h4666_dr60_u015", "PRJ-4666", "eo_0_633_sigv_40_CSR_0_150_Tau_0_.csv",
         40.0, 0.633, 0.15, 0.000, 20.0, "double_amplitude", 12.6953, 4),
    Case("h4666_dr90_u035", "PRJ-4666", "eo_0_545_sigv_40_CSR_0_350_Tau_0_.csv",
         40.0, 0.545, 0.35, 0.000, 28.0, "double_amplitude", 18.2188, 7),
    Case("h4666_dr50_b015_c018", "PRJ-4666", "eo_0_655_sigv_40_CSR_0_180_Tau_6_.csv",
         40.0, 0.655, 0.18, 0.150, 31.0, "single_amplitude", 20.1094, 8),
    Case("h4666_dr50_b015_c021", "PRJ-4666", "eo_0_657_sigv_40_CSR_0_210_Tau_6_.csv",
         40.0, 0.657, 0.21, 0.150, 16.0, "single_amplitude", 10.1172, 3),
    Case("h4666_dr60_b015_c022", "PRJ-4666", "eo_0_634_sigv_40_CSR_0_220_Tau_6_.csv",
         40.0, 0.634, 0.22, 0.150, 19.0, "single_amplitude", 12.1328, 4),
    Case("h4666_dr50_b0075_c017", "PRJ-4666", "eo_0_657_sigv_40_CSR_0_170_Tau_3_.csv",
         40.0, 0.657, 0.17, 0.075, 16.0, "single_amplitude", 10.1250, 3),
    Case("h4666_dr50_b025_c030", "PRJ-4666", "eo_0_660_sigv_40_CSR_0_300_Tau_10_.csv",
         40.0, 0.660, 0.30, 0.250, 14.0, "single_amplitude", 9.05469, 3),
    Case("h4666_dr60_b025_c033", "PRJ-4666", "eo_0_634_sigv_40_CSR_0_330_Tau_10_.csv",
         40.0, 0.634, 0.33, 0.250, 11.0, "single_amplitude", 7.14063, 2),
    Case("h4666_dr90_b025_c045", "PRJ-4666", "eo_0_553_sigv_40_CSR_0_450_Tau_10_.csv",
         40.0, 0.553, 0.45, 0.250, 35.0, "single_amplitude", 23.2422, 9),
    Case("h4666_dr50_b0375_c035", "PRJ-4666", "eo_0_662_sigv_40_CSR_0_350_Tau_15_.csv",
         40.0, 0.662, 0.35, 0.375, 20.0, "single_amplitude", 13.1875, 5),
    Case("h4666_dr90_b0375_c050", "PRJ-4666", "eo_0_550_sigv_40_CSR_0_500_Tau_15_.csv",
         40.0, 0.550, 0.50, 0.375, 44.0, "single_amplitude", 29.1875, 11),
)


# The filenames happen to be unique within the selected project/dataset pairs.
# Resolve them explicitly and reject ambiguity so later dataset additions cannot
# silently change a validation source.
def holdout_experiment_path(case: Case) -> Path:
    matches = sorted((WORKSPACE / case.dataset).rglob(case.filename))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one holdout file for {case.case_id}; found {len(matches)}"
        )
    return matches[0]


def validate_holdout_manifest() -> None:
    """Reject a changed or incorrectly transcribed experimental manifest."""
    for case in HOLDOUT_CASES:
        data = strict.load_designsafe_csv(holdout_experiment_path(case))
        scalar_checks = (
            ("void_ratio", case.void_ratio),
            ("vertical_effective_stress_0", case.vertical_stress),
            ("CSR", case.csr),
            ("static_shear_ratio", case.bias),
        )
        for name, expected in scalar_checks:
            if not np.isclose(float(data[name]), expected, atol=5.0e-4):
                raise ValueError(
                    f"{case.case_id}: {name}={data[name]} does not match {expected}"
                )
        gamma = np.asarray(data["gamma_percent"], dtype=float)
        gamma -= gamma[0]
        if case.criterion == "single_amplitude":
            measure = np.maximum.accumulate(np.abs(gamma))
        else:
            measure = (
                np.maximum.accumulate(gamma) - np.minimum.accumulate(gamma)
            )
        crossing = np.flatnonzero(measure >= 7.5)
        if crossing.size == 0:
            raise ValueError(f"{case.case_id}: experimental strain criterion is absent")
        measured_cycle = float(np.asarray(data["cycle"])[crossing[0]])
        if not np.isclose(measured_cycle, case.experimental_cycles, atol=6.0e-4):
            raise ValueError(
                f"{case.case_id}: N={measured_cycle} does not match "
                f"manifest N={case.experimental_cycles}"
            )


MODEL_TYPES = {
    "production": (RIVASandModel, RIVASandParameters),
    "signed_pt": (
        RIVASandLooseUnbiasedCorrectionModel,
        RIVASandLooseUnbiasedCorrectionParameters,
    ),
    "successor": (
        RIVASandLooseBiasedShearFlowModel,
        RIVASandLooseBiasedShearFlowParameters,
    ),
}


def configure_model(args: argparse.Namespace) -> None:
    model_type, parameter_type = MODEL_TYPES[args.model]
    strict.RIVASandDBFBModel = model_type
    strict.RIVASandDBFBParameters = parameter_type

    def parameters_from_args(_: argparse.Namespace):
        defaults = (
            asdict(reference_parameters())
            if args.model == "production"
            else asdict(parameter_type())
        )
        allowed = {item.name for item in fields(parameter_type)}
        overrides = strict._parameter_overrides(args.parameter, defaults)
        unknown = sorted(set(overrides) - allowed)
        if unknown:
            raise ValueError(
                f"unknown {args.model} parameter(s): {', '.join(unknown)}"
            )
        defaults.update(overrides)
        return parameter_type(**defaults)

    strict.parameters_from_args = parameters_from_args


def build_parser() -> argparse.ArgumentParser:
    parser = strict.build_parser()
    parser.description = __doc__
    parser.add_argument(
        "--model", choices=tuple(MODEL_TYPES), default="successor",
        help="constitutive checkpoint to evaluate",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    strict.CASES = HOLDOUT_CASES
    strict.experiment_path = holdout_experiment_path
    validate_holdout_manifest()
    configure_model(args)
    result = strict.audit(args)
    print(strict.json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
