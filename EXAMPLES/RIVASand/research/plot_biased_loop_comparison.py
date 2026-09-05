"""Plot production RIVA-Sand and the biased-loop prototype against DSS data."""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

import RIVASandBaselineAudit as audit  # noqa: E402


RESULTS = HERE / "results"
BASELINE = RESULTS / "baseline_full"
RESEARCH = RESULTS / "candidate_c1_full"
OUTPUT = RESULTS / "biased_loop_comparison"

EXPERIMENT_COLOR = "#D55E00"
BASELINE_COLOR = "#5B5B5B"
RESEARCH_COLOR = "#0072B2"
CASE_IDS = ("3484_b025", "4666_dense_b025", "4666_dense_b0375")


def read_history(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {
        name: np.atleast_1d(np.asarray(data[name], dtype=float))
        for name in data.dtype.names or ()
    }


def read_metrics(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return {row["case_id"]: row for row in csv.DictReader(stream)}


def cycle_mask(history: dict[str, np.ndarray], cycle: int) -> np.ndarray:
    return (history["cycle"] >= cycle) & (history["cycle"] < cycle + 1)


def plot_histories() -> None:
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    baseline_metrics = read_metrics(BASELINE / "metrics.csv")
    research_metrics = read_metrics(RESEARCH / "metrics.csv")
    fig, axes = plt.subplots(3, 3, figsize=(13.2, 10.4))
    for row, case_id in enumerate(CASE_IDS):
        case = cases[case_id]
        experiment = audit.audit_module.shifted_experiment(case)
        baseline = read_history(BASELINE / f"{case_id}.csv")
        research = read_history(RESEARCH / f"{case_id}.csv")
        end = min(12.0, float(experiment["cycle"][-1]))
        histories = (
            (experiment, "Experiment", EXPERIMENT_COLOR, "-", 1.35),
            (baseline, "Production", BASELINE_COLOR, ":", 1.20),
            (research, "Biased-loop research", RESEARCH_COLOR, "--", 1.20),
        )
        for history, label, color, style, width in histories:
            mask = history["cycle"] <= min(end, float(history["cycle"][-1]))
            ru_key = "ru" if label == "Experiment" else "ru_vertical"
            axes[row, 0].plot(
                history["gamma_percent"][mask], history["tau"][mask],
                color=color, ls=style, lw=width, label=label,
            )
            axes[row, 1].plot(
                history["vertical_effective_stress"][mask], history["tau"][mask],
                color=color, ls=style, lw=width,
            )
            axes[row, 2].plot(
                history["cycle"][mask], history[ru_key][mask],
                color=color, ls=style, lw=width,
            )
        base = baseline_metrics[case_id]
        new = research_metrics[case_id]
        axes[row, 0].set_ylabel(f"{case_id}\nShear stress, τ (kPa)")
        axes[row, 1].set_title(
            rf"$\alpha$={case.bias:g}, CSR={case.csr:g}; "
            f"range ratio {float(base['range_ratio']):.2f}→"
            f"{float(new['range_ratio']):.2f}"
        )
        for axis in axes[row]:
            axis.grid(color="0.88", lw=0.6)
    axes[0, 0].set_title("Stress–strain histories")
    axes[0, 2].set_title(r"Pore-pressure ratio, $r_u$")
    axes[-1, 0].set_xlabel(r"Shear strain, $\gamma$ (%)")
    axes[-1, 1].set_xlabel(r"Vertical effective stress, $\sigma'_v$ (kPa)")
    axes[-1, 2].set_xlabel("Cycle, N")
    axes[0, 0].legend(frameon=False, ncol=1, loc="best")
    fig.suptitle(
        "Biased DSS comparison: production RIVA-Sand and compact loop prototype",
        fontsize=14,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.975))
    fig.savefig(OUTPUT / "biased_history_comparison.png", dpi=240)
    fig.savefig(OUTPUT / "biased_history_comparison.pdf")
    plt.close(fig)


def plot_isolated_loops() -> None:
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1))
    for axis, case_id in zip(axes, CASE_IDS):
        case = cases[case_id]
        histories = (
            (audit.audit_module.shifted_experiment(case), "Experiment", EXPERIMENT_COLOR, "-"),
            (read_history(BASELINE / f"{case_id}.csv"), "Production", BASELINE_COLOR, ":"),
            (read_history(RESEARCH / f"{case_id}.csv"), "Biased-loop research", RESEARCH_COLOR, "--"),
        )
        for history, label, color, style in histories:
            mask = cycle_mask(history, case.comparison_cycle)
            axis.plot(
                history["gamma_percent"][mask], history["tau"][mask],
                color=color, ls=style, lw=1.7, label=label,
            )
        axis.set(
            title=(
                f"{case_id}: cycle {case.comparison_cycle}\n"
                rf"$\alpha$={case.bias:g}, CSR={case.csr:g}"
            ),
            xlabel=r"Shear strain, $\gamma$ (%)",
            ylabel=r"Shear stress, $\tau$ (kPa)",
        )
        axis.grid(color="0.88", lw=0.6)
    axes[0].legend(frameon=False, loc="best")
    fig.suptitle("Fixed-cycle biased stress–strain loop comparison", fontsize=14)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    fig.savefig(OUTPUT / "biased_isolated_loop_comparison.png", dpi=240)
    fig.savefig(OUTPUT / "biased_isolated_loop_comparison.pdf")
    plt.close(fig)


def write_summary() -> None:
    baseline = read_metrics(BASELINE / "metrics.csv")
    research = read_metrics(RESEARCH / "metrics.csv")
    fields = (
        "case_id", "quantity", "experiment", "production", "research"
    )
    rows = []
    mappings = (
        ("loop_center_percent", "center_experiment", "center_model"),
        ("loop_range_percent", "range_experiment", "range_model"),
        ("loop_area_kPa_percent", "area_experiment", "area_model"),
        ("damping_ratio", "damping_experiment", "damping_model"),
        ("cycles_to_criterion", "N_experiment", "N_model"),
    )
    for case_id in CASE_IDS:
        for quantity, experiment_key, model_key in mappings:
            rows.append({
                "case_id": case_id,
                "quantity": quantity,
                "experiment": baseline[case_id][experiment_key],
                "production": baseline[case_id][model_key],
                "research": research[case_id][model_key],
            })
    with (OUTPUT / "biased_loop_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_histories()
    plot_isolated_loops()
    write_summary()
    print(OUTPUT)


if __name__ == "__main__":
    main()
