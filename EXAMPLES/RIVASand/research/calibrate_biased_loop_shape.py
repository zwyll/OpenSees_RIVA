"""Bounded loop-shape calibration for the compact biased-loop prototype.

The search is intentionally small.  It does not calibrate pore pressure or
cycles to liquefaction and it never changes the production RIVA-Sand oracle.
Candidates are screened at 16 points/cycle with two constitutive substeps;
only the four best are reevaluated at 32 points/cycle and four substeps.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
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
from RIVASandBiasedLoop import (  # noqa: E402
    RIVASandBiasedLoopModel,
    RIVASandBiasedLoopParameters,
)


baseline.audit_module.RIVASandDBFBModel = RIVASandBiasedLoopModel
baseline.audit_module.RIVASandDBFBParameters = RIVASandBiasedLoopParameters


CASES = tuple(
    case for case in baseline.audit_module.CASES
    if case.case_id in {"3484_b025", "4666_dense_b025", "4666_dense_b0375"}
)
OUTPUT = HERE / "results" / "biased_loop_calibration"


def score_record(record: dict[str, Any]) -> float:
    if not int(record["comparison_cycle_available"]):
        return 25.0
    normalized = (
        (0.20, float(record["raw_phase_nrmse"]) / 0.25),
        (0.20, float(record["affine_shape_rmse"]) / 0.18),
        (0.12, math.log(max(float(record["range_ratio"]), 1.0e-8)) / math.log(1.30)),
        (0.10, math.log(max(float(record["area_ratio"]), 1.0e-8)) / math.log(1.50)),
        (0.10, float(record["center_normalized_error"]) / 0.25),
        (0.06, float(record["damping_absolute_error"]) / 0.05),
        (0.10, float(record["mid_opening_error"]) / 0.15),
        (0.08, float(record["center_trajectory_nrmse"]) / 0.30),
        (
            0.04,
            max(float(record["negative_center_drift_ratio"]) - 0.10, 0.0)
            / 0.10,
        ),
    )
    return float(math.sqrt(sum(weight * value * value for weight, value in normalized)))


def evaluate(task: tuple[int, dict[str, float], int, int]) -> dict[str, Any]:
    index, values, points, substeps = task
    parameters = replace(RIVASandBiasedLoopParameters(), **values)
    records = []
    for case in CASES:
        record, _ = baseline.audit_module.compare_case(
            parameters, case, points=points, substeps=substeps, full=False
        )
        records.append(record)
    return {
        "index": index,
        "parameters": values,
        "points_per_cycle": points,
        "substeps": substeps,
        "score": float(np.mean([score_record(record) for record in records])),
        "records": records,
    }


def candidates() -> list[dict[str, float]]:
    result = [
        {},
        {
            "cyclic_shear_modulus_reduction": 0.87,
            "branch_compliance_gain": 4.5,
            "branch_compliance_exponent": 1.0,
            "branch_directional_balance": 0.03,
        },
        {
            "cyclic_shear_modulus_reduction": 0.85,
            "branch_compliance_gain": 5.0,
            "branch_compliance_exponent": 1.0,
            "branch_directional_balance": 0.03,
        },
    ]
    rng = np.random.default_rng(20260827)
    for _ in range(29):
        result.append({
            "cyclic_shear_modulus_reduction": float(rng.uniform(0.84, 0.90)),
            "branch_compliance_gain": float(rng.uniform(1.5, 7.0)),
            "branch_compliance_exponent": float(rng.uniform(0.45, 1.8)),
            "branch_directional_balance": float(rng.uniform(-0.18, 0.08)),
            "branch_balance_bias_exponent": float(rng.uniform(0.0, 0.8)),
            "branch_compliance_bias_exponent": float(rng.uniform(0.5, 1.8)),
        })
    return result


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    values = candidates()
    low_tasks = [(index, value, 16, 2) for index, value in enumerate(values)]
    with ProcessPoolExecutor(max_workers=4) as executor:
        low = list(executor.map(evaluate, low_tasks))
    low.sort(key=lambda item: item["score"])
    high_tasks = [
        (item["index"], item["parameters"], 32, 4) for item in low[:4]
    ]
    with ProcessPoolExecutor(max_workers=4) as executor:
        high = list(executor.map(evaluate, high_tasks))
    high.sort(key=lambda item: item["score"])
    payload = {
        "scope": "stress-strain loop shape only; Ru and CSR-N excluded",
        "default_parameters": asdict(RIVASandBiasedLoopParameters()),
        "low_resolution_ranking": low,
        "high_resolution_ranking": high,
    }
    path = OUTPUT / "bounded_search.json"
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n")
    print(json.dumps({
        "best_score": high[0]["score"],
        "best_parameters": high[0]["parameters"],
        "output": str(path),
    }, indent=2))


if __name__ == "__main__":
    main()
