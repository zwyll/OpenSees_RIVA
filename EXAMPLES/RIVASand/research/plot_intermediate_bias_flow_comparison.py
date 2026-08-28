"""Plot the two PRJ-3484 loops targeted by intermediate biased flow."""

from __future__ import annotations

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

import RIVASandBaselineAudit as baseline  # noqa: E402


strict = baseline.audit_module
RESULTS = HERE / "results" / "intermediate_bias_flow"
OLD = WORKSPACE / "rivasand_model_comparison" / "results_32x4" / "histories" / "current"
SHEAR = {
    case_id: RESULTS / "final_original2_32x4" / f"{case_id}.csv"
    for case_id in ("3484_b025", "3484_b030")
}
NEW = {
    case_id: RESULTS / "volumetric_final_target_32x4" / f"{case_id}.csv"
    for case_id in ("3484_b025", "3484_b030")
}
COLORS = {
    "experiment": "#D55E00",
    "old": "#777777",
    "shear": "#7B3294",
    "new": "#0072B2",
}


def load(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {
        name: np.atleast_1d(np.asarray(data[name], dtype=float))
        for name in data.dtype.names or ()
    }


def complete_cycle(history: dict[str, np.ndarray], cycle: int) -> np.ndarray:
    return (
        (history["cycle"] >= cycle - 1.0e-10)
        & (history["cycle"] <= cycle + 1.0 + 1.0e-10)
    )


def centers(history: dict[str, np.ndarray], maximum: int) -> tuple[np.ndarray, np.ndarray]:
    cycles: list[float] = []
    values: list[float] = []
    for cycle in range(maximum + 1):
        mask = complete_cycle(history, cycle)
        if np.count_nonzero(mask) < 8:
            continue
        gamma = history["gamma_percent"][mask]
        values.append(0.5 * (float(np.min(gamma)) + float(np.max(gamma))))
        cycles.append(cycle + 0.5)
    return np.asarray(cycles), np.asarray(values)


def main() -> None:
    cases = [
        next(item for item in strict.CASES if item.case_id == case_id)
        for case_id in ("3484_b025", "3484_b030")
    ]
    fig, axes = plt.subplots(2, 4, figsize=(16.0, 7.8), squeeze=False)
    for row, case in enumerate(cases):
        experiment = strict.shifted_experiment(case)
        old = load(OLD / f"{case.case_id}.csv")
        shear = load(SHEAR[case.case_id])
        new = load(NEW[case.case_id])
        cycle = case.comparison_cycle
        histories = (
            (experiment, "experiment", "Experiment", "-", "ru"),
            (old, "old", "Previous mapping/backstress", "--", "ru_vertical"),
            (shear, "shear", "Shear-flow correction", "-.", "ru_vertical"),
            (new, "new", "Coupled volume correction", ":", "ru_vertical"),
        )
        for history, key, label, style, ru_name in histories:
            mask = complete_cycle(history, cycle)
            axes[row, 0].plot(
                history["gamma_percent"][mask], history["tau"][mask],
                color=COLORS[key], ls=style, lw=2.0, label=label,
            )
            axes[row, 1].plot(
                history["vertical_effective_stress"][mask],
                history["tau"][mask], color=COLORS[key], ls=style, lw=2.0,
            )
            end = min(cycle + 1.0, float(history["cycle"][-1]))
            history_mask = history["cycle"] <= end + 1.0e-10
            axes[row, 2].plot(
                history["cycle"][history_mask], history[ru_name][history_mask],
                color=COLORS[key], ls=style, lw=1.8,
            )
            center_cycle, center = centers(history, cycle)
            axes[row, 3].plot(
                center_cycle, center, color=COLORS[key], ls=style, lw=1.8,
                marker="o", ms=3.0,
            )
        axes[row, 0].set_ylabel(
            f"{case.case_id}, cycle {cycle}\n" + r"$\tau$ (kPa)"
        )
        axes[row, 0].set_xlabel(r"Shear strain, $\gamma$ (%)")
        axes[row, 1].set_xlabel(r"Vertical effective stress, $\sigma'_v$ (kPa)")
        axes[row, 1].set_ylabel(r"$\tau$ (kPa)")
        axes[row, 2].set_xlabel("Cycle, N")
        axes[row, 2].set_ylabel(r"$r_u$")
        axes[row, 3].set_xlabel("Cycle, N")
        axes[row, 3].set_ylabel(r"Loop center, $\gamma_c$ (%)")
        axes[row, 0].text(
            0.03, 0.95,
            rf"$N_{{exp}}={case.experimental_cycles:.2f}$" + "\n"
            + rf"$N_{{new}}={float(new['cycle'][-1]):.2f}$",
            transform=axes[row, 0].transAxes, ha="left", va="top", fontsize=9,
        )
        for axis in axes[row]:
            axis.grid(color="#D9D9D9", lw=0.65)
            axis.set_axisbelow(True)
    titles = (
        "Complete cycle-matched stress–strain loop",
        "Cycle-matched effective-stress path",
        "Pore-pressure history",
        "Accumulated loop-center drift",
    )
    for axis, title in zip(axes[0], titles):
        axis.set_title(title)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=3, frameon=False,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "PRJ-3484 intermediate-density biased tests: shear and volume corrections",
        y=1.015, fontsize=15,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
    RESULTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(RESULTS / "intermediate_bias_flow_before_after.png", dpi=260, bbox_inches="tight")
    fig.savefig(RESULTS / "intermediate_bias_flow_before_after.pdf", bbox_inches="tight")

    old_csr = load(RESULTS / "csr_final_32x4" / "csr_records.csv")
    new_csr = load(RESULTS / "volumetric_final_csr_32x4" / "csr_records.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), squeeze=False)
    for axis, stress, alpha in zip(axes[0], (40.0, 100.0), (0.25, 0.30)):
        old_mask = np.isclose(old_csr["vertical_stress"], stress)
        new_mask = np.isclose(new_csr["vertical_stress"], stress)
        order_old = np.argsort(old_csr["CSR"][old_mask])
        order_new = np.argsort(new_csr["CSR"][new_mask])
        x_old = old_csr["CSR"][old_mask][order_old]
        x_new = new_csr["CSR"][new_mask][order_new]
        duration = new_csr["experiment_duration"][new_mask][order_new]
        exp_n = new_csr["N_experiment"][new_mask][order_new]
        old_n = old_csr["N_model"][old_mask][order_old]
        new_n = new_csr["N_model"][new_mask][order_new]
        for values, color, marker, label in (
            (exp_n, COLORS["experiment"], "o", "Experiment"),
            (old_n, COLORS["shear"], "^", "Shear-flow correction"),
            (new_n, COLORS["new"], "s", "Coupled volume correction"),
        ):
            x = x_new if label != "Shear-flow correction" else x_old
            runout_duration = duration
            plotted = np.where(np.isfinite(values), values, runout_duration)
            finite = np.isfinite(values)
            axis.plot(x, plotted, color=color, lw=1.8, ls="--" if label != "Experiment" else "-")
            axis.scatter(x[finite], plotted[finite], color=color, marker=marker, s=42, label=label)
            axis.scatter(
                x[~finite], plotted[~finite], facecolors="white", edgecolors=color,
                marker=marker, s=48, linewidths=1.6,
            )
        axis.set_yscale("log")
        axis.set_xlabel("CSR")
        axis.set_ylabel(r"Cycles to 7.5% DA, $N_{7.5\%DA}$")
        axis.set_title(rf"$\sigma'_{{v0}}={stress:.0f}$ kPa, $\alpha={alpha:.2f}$")
        axis.grid(color="#D9D9D9", lw=0.65, which="both")
        axis.set_axisbelow(True)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=3, frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.suptitle("PRJ-3484 CSR–N audit of the affected bias windows", y=0.91)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.80))
    fig.savefig(RESULTS / "intermediate_bias_flow_csr_n.png", dpi=260, bbox_inches="tight")
    fig.savefig(RESULTS / "intermediate_bias_flow_csr_n.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
