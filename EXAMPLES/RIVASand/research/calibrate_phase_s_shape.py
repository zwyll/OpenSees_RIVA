"""High-resolution bounded search for biased DSS S-shaped loop curvature."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

import RIVASandBaselineAudit as baseline  # noqa: E402
from RIVASandPhaseTransformation import (  # noqa: E402
    RIVASandPhaseTransformationModel,
    RIVASandPhaseTransformationParameters,
)


baseline.audit_module.RIVASandDBFBModel = RIVASandPhaseTransformationModel
baseline.audit_module.RIVASandDBFBParameters = RIVASandPhaseTransformationParameters
CASES = tuple(
    case for case in baseline.audit_module.CASES
    if case.case_id in {"4666_dense_b025", "4666_dense_b0375"}
)
OUTPUT = HERE / "results" / "phase_s_shape_calibration"


def score(record: dict[str, Any]) -> float:
    if not int(record["comparison_cycle_available"]):
        return 30.0
    terms = (
        (0.24, float(record["affine_shape_rmse"]) / 0.18),
        (0.16, float(record["raw_phase_nrmse"]) / 0.25),
        (0.12, float(record["mid_opening_error"]) / 0.15),
        (0.10, math.log(max(float(record["range_ratio"]), 1e-8)) / math.log(1.30)),
        (0.08, math.log(max(float(record["area_ratio"]), 1e-8)) / math.log(1.50)),
        (0.08, float(record["damping_absolute_error"]) / 0.05),
        (0.08, float(record["center_normalized_error"]) / 0.25),
        (0.06, float(record["center_trajectory_nrmse"]) / 0.30),
        (0.05, float(record["ru_mean_error"]) / 0.10),
        (0.03, float(record["ru_rmse"]) / 0.12),
    )
    return float(math.sqrt(sum(weight * value * value for weight, value in terms)))


def evaluate(task: tuple[int, dict[str, float]]) -> dict[str, Any]:
    index, values = task
    parameters = replace(RIVASandPhaseTransformationParameters(), **values)
    records = []
    for case in CASES:
        record, _ = baseline.audit_module.compare_case(
            parameters, case, points=32, substeps=4, full=False
        )
        records.append(record)
    return {
        "index": index,
        "parameters": values,
        "score": float(np.mean([score(record) for record in records])),
        "records": records,
    }


def candidates() -> list[dict[str, float]]:
    values = [
        {},
        {
            "cyclic_shear_modulus_reduction": 0.78,
            "phase_compliance_peak": 3.0,
            "phase_compliance_bell_gain": 5.0,
        },
        {
            "cyclic_shear_modulus_reduction": 0.72,
            "phase_compliance_peak": 4.0,
            "phase_compliance_bell_gain": 7.0,
        },
    ]
    rng = np.random.default_rng(20260828)
    for _ in range(33):
        values.append({
            "cyclic_shear_modulus_reduction": float(rng.uniform(0.60, 0.84)),
            "phase_compliance_peak": float(rng.uniform(1.5, 6.0)),
            "phase_compliance_bell_gain": float(rng.uniform(1.0, 10.0)),
            "phase_compliance_shape": float(rng.uniform(0.5, 1.6)),
            "phase_compliance_location": float(rng.uniform(0.35, 0.55)),
            "phase_compliance_half_width": float(rng.uniform(0.50, 0.85)),
            "branch_directional_balance": float(rng.uniform(0.025, 0.055)),
            "branch_balance_bias_exponent": float(rng.uniform(0.4, 1.0)),
        })
    return values


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    tasks = list(enumerate(candidates()))
    with ProcessPoolExecutor(max_workers=4) as executor:
        ranking = list(executor.map(evaluate, tasks))
    ranking.sort(key=lambda item: item["score"])
    path = OUTPUT / "bounded_search.json"
    path.write_text(json.dumps({
        "scope": "two dense biased histories at 32 points/cycle and 4 substeps",
        "ranking": ranking,
    }, indent=2, allow_nan=True) + "\n")
    print(json.dumps({
        "best_score": ranking[0]["score"],
        "best_parameters": ranking[0]["parameters"],
        "output": str(path),
    }, indent=2))


if __name__ == "__main__":
    main()
