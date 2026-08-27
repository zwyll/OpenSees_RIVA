"""Plot the late-cycle accumulation-control checkpoint against DSS data."""

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
PRIOR = RESULTS / "phase_transformation_density_extended_final_full"
CURRENT = RESULTS / "accumulation_control_host_committed_final_full"
OUTPUT = RESULTS / "accumulation_control_comparison"

EXPERIMENT_COLOR = "#D55E00"
PRIOR_COLOR = "#777777"
CURRENT_COLOR = "#0072B2"


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


def double_amplitude(history: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    cycle = np.concatenate(([0.0], history["cycle"]))
    strain = np.concatenate(([0.0], history["gamma_percent"]))
    return cycle, np.maximum.accumulate(strain) - np.minimum.accumulate(strain)


def loop_center(history: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    last_complete = int(np.floor(float(history["cycle"][-1])))
    cycles, centers = [], []
    for cycle in range(1, last_complete + 1):
        mask = cycle_mask(history, cycle)
        if np.count_nonzero(mask) < 3:
            continue
        gamma = history["gamma_percent"][mask]
        cycles.append(float(cycle))
        centers.append(0.5 * (float(np.max(gamma)) + float(np.min(gamma))))
    return np.asarray(cycles), np.asarray(centers)


def case_histories(case_id: str):
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    case = cases[case_id]
    return case, (
        (audit.audit_module.shifted_experiment(case), "Experiment", EXPERIMENT_COLOR, "-", 2.0),
        (read_history(PRIOR / f"{case_id}.csv"), "Prior PT", PRIOR_COLOR, ":", 1.7),
        (read_history(CURRENT / f"{case_id}.csv"), "Restrained accumulation", CURRENT_COLOR, "--", 1.9),
    )


def plot_timing_and_early_loops() -> None:
    metrics_prior = read_metrics(PRIOR / "metrics.csv")
    metrics_current = read_metrics(CURRENT / "metrics.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4))
    for column, case_id in enumerate(("3484_b025", "3484_b030")):
        case, histories = case_histories(case_id)
        for history, label, color, style, width in histories:
            cycle, measure = double_amplitude(history)
            axes[0, column].plot(
                cycle, measure, color=color, ls=style, lw=width, label=label,
            )
        axes[0, column].axhline(7.5, color="black", lw=1.0, alpha=0.7)
        n_values = (
            case.experimental_cycles,
            float(metrics_prior[case_id]["N_model"]),
            float(metrics_current[case_id]["N_model"]),
        )
        axes[0, column].scatter(
            n_values, [7.5] * 3,
            color=(EXPERIMENT_COLOR, PRIOR_COLOR, CURRENT_COLOR),
            s=34, zorder=5,
        )
        axes[0, column].set(
            title=(
                f"{case_id}: 7.5% DA criterion\n"
                f"Experiment {n_values[0]:.3g}, prior {n_values[1]:.3g}, "
                f"restrained {n_values[2]:.3g} cycles"
            ),
            xlabel="Cycle, N",
            ylabel="Double-amplitude shear strain (%)",
            ylim=(0.0, 8.2),
        )
        axes[0, column].set_xlim(
            0.0, 70.0 if case_id == "3484_b025" else 3.5
        )

        for history, label, color, style, width in histories:
            mask = cycle_mask(history, case.comparison_cycle)
            axes[1, column].plot(
                history["gamma_percent"][mask], history["tau"][mask],
                color=color, ls=style, lw=width, label=label,
            )
        axes[1, column].set(
            title=f"{case_id}: cycle-{case.comparison_cycle} early loop",
            xlabel=r"Shear strain, $\gamma$ (%)",
            ylabel=r"Shear stress, $\tau$ (kPa)",
        )

    axes[0, 0].legend(frameon=False, fontsize=9, loc="upper left")
    for axis in axes.flat:
        axis.grid(color="0.88", lw=0.6)
    fig.suptitle(
        "Late-cycle restraint delays ratcheting without changing the calibrated early loop"
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(OUTPUT / "accumulation_timing_and_early_loops.png", dpi=300)
    fig.savefig(OUTPUT / "accumulation_timing_and_early_loops.pdf")
    plt.close(fig)


def plot_b025_evolution() -> None:
    case, histories = case_histories("3484_b025")
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.7))
    for history, label, color, style, width in histories:
        cycles, centers = loop_center(history)
        axes[0].plot(cycles, centers, color=color, ls=style, lw=width, label=label)
        ru_key = "ru" if label == "Experiment" else "ru_vertical"
        axes[1].plot(
            history["cycle"], history[ru_key], color=color, ls=style, lw=width,
        )
        mask = cycle_mask(history, case.comparison_cycle)
        axes[2].plot(
            history["vertical_effective_stress"][mask], history["tau"][mask],
            color=color, ls=style, lw=width,
        )
    axes[0].set(
        title="Cycle-by-cycle strain-center drift",
        xlabel="Cycle, N", ylabel=r"Loop center, $\gamma_c$ (%)", xlim=(0.0, 70.0),
    )
    axes[1].set(
        title=r"Pore-pressure ratio history",
        xlabel="Cycle, N", ylabel=r"$r_u$", xlim=(0.0, 70.0),
    )
    axes[2].set(
        title="Cycle-11 effective-stress path",
        xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
        ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[0].legend(frameon=False, fontsize=9)
    for axis in axes:
        axis.grid(color="0.88", lw=0.6)
    fig.suptitle("PRJ-3484 b025: isolated effect of the late-cycle shakedown state")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(OUTPUT / "3484_b025_accumulation_evolution.png", dpi=300)
    fig.savefig(OUTPUT / "3484_b025_accumulation_evolution.pdf")
    plt.close(fig)


def plot_objectivity() -> None:
    coarse = read_metrics(CURRENT / "metrics.csv")
    refined = read_metrics(
        RESULTS / "accumulation_control_host_committed_64x2_four" / "metrics.csv"
    )
    case_ids = ("3484_b025", "3484_b030", "4666_dense_b025", "4666_dense_b0375")
    x = np.arange(len(case_ids), dtype=float)
    width = 0.34
    n_coarse = [float(coarse[key]["N_model"]) for key in case_ids]
    n_refined = [float(refined[key]["N_model"]) for key in case_ids]
    fig, axis = plt.subplots(figsize=(9.0, 4.6))
    axis.bar(x - width / 2.0, n_coarse, width, color=CURRENT_COLOR, label="32 points/cycle, 4 substeps")
    axis.bar(x + width / 2.0, n_refined, width, color="#E69F00", label="64 points/cycle, 2 substeps")
    for index, (left, right) in enumerate(zip(n_coarse, n_refined)):
        difference = 100.0 * (right / left - 1.0)
        axis.text(index, max(left, right) + 0.8, f"{difference:+.1f}%", ha="center", fontsize=9)
    axis.set(
        title="DSS timestep/substep audit: cycles to 7.5% DA strain",
        ylabel="Cycles, N", xticks=x, xticklabels=case_ids,
    )
    axis.grid(axis="y", color="0.88", lw=0.6)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT / "accumulation_objectivity.png", dpi=300)
    fig.savefig(OUTPUT / "accumulation_objectivity.pdf")
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_timing_and_early_loops()
    plot_b025_evolution()
    plot_objectivity()
    print(OUTPUT)


if __name__ == "__main__":
    main()
