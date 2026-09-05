"""Bounded joint calibration of biased loop shape and pore-pressure response.

The search is intentionally limited to the three biased calibration histories.
It changes neither the production source nor the frozen validation histories.
Candidates that fail to reach a declared comparison cycle receive a dominant
penalty; liquefaction cycles and all five non-calibration cases are evaluated
only after the bounded search.
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
from RIVASandPhaseTransformation import (  # noqa: E402
    RIVASandPhaseTransformationModel,
    RIVASandPhaseTransformationParameters,
)


baseline.audit_module.RIVASandDBFBModel = RIVASandPhaseTransformationModel
baseline.audit_module.RIVASandDBFBParameters = RIVASandPhaseTransformationParameters

CASES = tuple(
    case for case in baseline.audit_module.CASES
    if case.case_id in {"3484_b025", "4666_dense_b025", "4666_dense_b0375"}
)
OUTPUT = HERE / "results" / "phase_transformation_calibration"


def score_record(record: dict[str, Any]) -> float:
    if not int(record["comparison_cycle_available"]):
        return 30.0
    normalized = (
        (0.12, float(record["raw_phase_nrmse"]) / 0.25),
        (0.10, float(record["affine_shape_rmse"]) / 0.18),
        (0.08, math.log(max(float(record["range_ratio"]), 1.0e-8)) / math.log(1.30)),
        (0.06, math.log(max(float(record["area_ratio"]), 1.0e-8)) / math.log(1.50)),
        (0.08, float(record["center_normalized_error"]) / 0.25),
        (0.05, float(record["damping_absolute_error"]) / 0.05),
        (0.07, float(record["mid_opening_error"]) / 0.15),
        (0.08, float(record["center_trajectory_nrmse"]) / 0.30),
        (0.16, float(record["ru_rmse"]) / 0.12),
        (0.10, float(record["ru_mean_error"]) / 0.10),
    )
    wave = 0.0
    if np.isfinite(float(record["ru_phase_error_deg"])):
        wave = math.hypot(
            math.log(max(float(record["ru_amplitude_ratio"]), 1.0e-8))
            / math.log(1.40),
            float(record["ru_phase_error_deg"]) / 30.0,
        )
    regularity = max(
        float(record["negative_center_drift_ratio"]) - 0.10, 0.0
    ) / 0.10
    total = sum(weight * value * value for weight, value in normalized)
    total += 0.06 * wave * wave + 0.04 * regularity * regularity
    return float(math.sqrt(total))


def evaluate(task: tuple[int, dict[str, float], int, int]) -> dict[str, Any]:
    index, values, points, substeps = task
    parameters = replace(RIVASandPhaseTransformationParameters(), **values)
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
            "phase_compliance_peak": 4.5,
            "phase_compliance_shape": 1.0,
            "phase_compliance_location": 0.50,
            "phase_compliance_half_width": 0.50,
            "phase_contraction_rate": 0.0005,
            "phase_reversible_scale": -0.0010,
        },
        {
            "phase_compliance_peak": 4.5,
            "phase_compliance_shape": 2.0,
            "phase_compliance_location": 0.45,
            "phase_compliance_half_width": 0.55,
            "phase_contraction_rate": 0.0004,
            "phase_reversible_scale": -0.0010,
        },
    ]
    rng = np.random.default_rng(20260827)
    for _ in range(61):
        result.append({
            "phase_compliance_peak": float(rng.uniform(2.5, 6.5)),
            "phase_compliance_shape": float(rng.uniform(0.65, 3.0)),
            "phase_compliance_location": float(rng.uniform(0.30, 0.65)),
            "phase_compliance_half_width": float(rng.uniform(0.35, 0.75)),
            "branch_directional_balance": float(rng.uniform(-0.04, 0.08)),
            "phase_ratio": float(rng.uniform(0.55, 0.85)),
            "phase_width": float(rng.uniform(0.12, 0.30)),
            "phase_contraction_rate": float(rng.uniform(0.00025, 0.00070)),
            "phase_reversible_scale": float(rng.uniform(-0.0014, -0.0006)),
            "phase_reversible_relaxation_strain": float(
                rng.uniform(0.0007, 0.0030)
            ),
            "phase_pressure_exponent": float(rng.uniform(1.0, 4.0)),
            "phase_bias_exponent": float(rng.uniform(0.5, 2.0)),
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
        (item["index"], item["parameters"], 32, 4) for item in low[:8]
    ]
    with ProcessPoolExecutor(max_workers=4) as executor:
        high = list(executor.map(evaluate, high_tasks))
    high.sort(key=lambda item: item["score"])
    payload = {
        "scope": "three biased histories; loop shape and Ru jointly",
        "default_parameters": asdict(RIVASandPhaseTransformationParameters()),
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
