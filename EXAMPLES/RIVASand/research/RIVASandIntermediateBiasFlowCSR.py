"""Affected CSR--N audit for the intermediate biased-flow successor."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import json
import math
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

import RIVASandBaselineAudit as baseline  # noqa: E402
from RIVASandIntermediateBiasFlow import (  # noqa: E402
    RIVASandIntermediateBiasFlowModel,
    RIVASandIntermediateBiasFlowParameters,
)


strict = baseline.audit_module
SOURCE = (
    WORKSPACE / "pj_model_v6" / "results" / "static_shear_development"
    / "v6_static_shear_all_tests.csv"
)


def jobs() -> list[dict[str, float | str]]:
    result: list[dict[str, float | str]] = []
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        for index, row in enumerate(csv.DictReader(stream)):
            stress = float(row["vertical_stress"])
            alpha = float(row["static_shear_ratio"])
            if not (
                (math.isclose(stress, 40.0) and math.isclose(alpha, 0.25))
                or (math.isclose(stress, 100.0) and math.isclose(alpha, 0.30))
            ):
                continue
            result.append({
                "job_id": f"csr_{index:02d}",
                "vertical_stress": stress,
                "static_shear_stress": float(row["static_shear_stress"]),
                "alpha": alpha,
                "CSR": float(row["CSR"]),
                "void_ratio": 0.601,
                "N_experiment": float(row["N_experiment"]),
                "experiment_duration": float(row["experiment_duration"]),
            })
    return result


def run_one(
    task: tuple[dict[str, float | str], int, int, dict[str, object]]
) -> tuple[dict[str, object], dict[str, object]]:
    values, points, substeps, parameter_values = task
    strict.RIVASandDBFBModel = RIVASandIntermediateBiasFlowModel
    strict.RIVASandDBFBParameters = RIVASandIntermediateBiasFlowParameters
    case = strict.Case(
        case_id=str(values["job_id"]),
        dataset="PRJ-3484",
        filename="",
        vertical_stress=float(values["vertical_stress"]),
        void_ratio=float(values["void_ratio"]),
        csr=float(values["CSR"]),
        bias=float(values["alpha"]),
        maximum_cycles=float(values["experiment_duration"]),
        criterion="double_amplitude",
        experimental_cycles=float(values["N_experiment"]),
        comparison_cycle=0,
    )
    record, history = strict.run_case(
        RIVASandIntermediateBiasFlowParameters(**parameter_values),
        case,
        points=points,
        substeps=substeps,
        full=True,
    )
    record.update(values)
    return record, history


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--points", type=int, default=32)
    parser.add_argument("--substeps", type=int, default=4)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--job", action="append", default=[], help="run one CSR job id"
    )
    parser.add_argument(
        "--parameter", action="append", default=[], metavar="NAME=VALUE",
        help="override an intermediate-flow parameter",
    )
    args = parser.parse_args()
    parameter_values = asdict(RIVASandIntermediateBiasFlowParameters())
    parameter_values.update(
        strict._parameter_overrides(args.parameter, parameter_values)
    )
    selected = set(args.job)
    selected_jobs = [
        item for item in jobs()
        if not selected or str(item["job_id"]) in selected
    ]
    unknown = selected - {str(item["job_id"]) for item in selected_jobs}
    if unknown:
        raise ValueError(f"unknown CSR job(s): {', '.join(sorted(unknown))}")
    tasks = [
        (item, args.points, args.substeps, parameter_values)
        for item in selected_jobs
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run_one, tasks))
    records = [item[0] for item in results]
    histories = {
        str(record["job_id"]): history for record, history in results
    }
    records.sort(key=lambda item: str(item["job_id"]))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "csr_records.json").write_text(
        json.dumps(records, allow_nan=True, indent=2) + "\n",
        encoding="utf-8",
    )
    with (args.output / "csr_records.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    for case_id, history in histories.items():
        if not history or not len(history.get("cycle", [])):
            continue
        with (args.output / f"{case_id}.csv").open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.writer(stream)
            writer.writerow(history)
            writer.writerows(zip(*(history[key] for key in history)))
    for record in records:
        print(
            record["job_id"],
            f"CSR={record['CSR']}",
            f"N_exp={record['N_experiment']}",
            f"N_model={record['N_model']}",
            record["status"],
            flush=True,
        )


if __name__ == "__main__":
    main()
