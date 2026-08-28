"""Regenerate the complete two-dataset DSS comparison for the latest model."""

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
LATEST = RESULTS / "loose_biased_shear_flow_final_full"
OUTPUT = RESULTS / "loose_biased_shear_flow_all_datasets_comparison"

EXPERIMENT_COLOR = "#D55E00"
NUMERICAL_COLOR = "#0072B2"
THRESHOLD_COLOR = "#222222"
DATASETS = ("PRJ-3484", "PRJ-4666")


def read_history(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {
        name: np.atleast_1d(np.asarray(data[name], dtype=float))
        for name in data.dtype.names or ()
    }


def read_metrics() -> dict[str, dict[str, str]]:
    with (LATEST / "metrics.csv").open(newline="", encoding="utf-8") as stream:
        return {row["case_id"]: row for row in csv.DictReader(stream)}


def cases_for(dataset: str):
    return tuple(case for case in audit.audit_module.CASES if case.dataset == dataset)


def histories(case):
    return (
        audit.audit_module.shifted_experiment(case),
        read_history(LATEST / f"{case.case_id}.csv"),
    )


def cycle_mask(history: dict[str, np.ndarray], cycle: int) -> np.ndarray:
    return (history["cycle"] >= cycle) & (history["cycle"] < cycle + 1)


def criterion_history(
    history: dict[str, np.ndarray], criterion: str
) -> tuple[np.ndarray, np.ndarray]:
    cycle = np.concatenate(([0.0], history["cycle"]))
    strain = np.concatenate(([0.0], history["gamma_percent"]))
    running_min = np.minimum.accumulate(strain)
    running_max = np.maximum.accumulate(strain)
    if criterion == "double_amplitude":
        measure = running_max - running_min
    elif criterion == "single_amplitude":
        measure = np.maximum(np.abs(running_min), np.abs(running_max))
    else:
        raise ValueError(f"unknown strain criterion {criterion}")
    return cycle, measure


def apply_grid(axis) -> None:
    axis.grid(color="0.88", lw=0.6)
    axis.set_axisbelow(True)


def dataset_slug(dataset: str) -> str:
    return dataset.lower().replace("-", "_")


def plot_history_matrix(dataset: str) -> None:
    cases = cases_for(dataset)
    fig, axes = plt.subplots(len(cases), 3, figsize=(14.0, 12.2), squeeze=False)
    for row, case in enumerate(cases):
        experiment, numerical = histories(case)
        for history, label, color, style, width, ru_key in (
            (experiment, "Experiment", EXPERIMENT_COLOR, "-", 1.7, "ru"),
            (numerical, "Corrected research successor", NUMERICAL_COLOR, "--", 1.5, "ru_vertical"),
        ):
            end = min(12.0, float(history["cycle"][-1]))
            mask = history["cycle"] <= end
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
        axes[row, 0].set_ylabel(
            f"{case.case_id}\n"
            rf"$\alpha={case.bias:g}$, CSR$={case.csr:g}$\n"
            r"$\tau$ (kPa)"
        )
        for axis in axes[row]:
            apply_grid(axis)
    axes[0, 0].set_title("Stress–strain history (first 12 cycles)")
    axes[0, 1].set_title("Measured DSS effective-stress path")
    axes[0, 2].set_title(r"Pore-pressure ratio, $r_u$")
    axes[-1, 0].set_xlabel(r"Shear strain, $\gamma$ (%)")
    axes[-1, 1].set_xlabel(r"Vertical effective stress, $\sigma'_v$ (kPa)")
    axes[-1, 2].set_xlabel("Cycle, N")
    axes[0, 0].legend(frameon=False, fontsize=9, ncol=2)
    fig.suptitle(f"{dataset}: corrected DSS research-successor comparison")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.975))
    stem = OUTPUT / f"{dataset_slug(dataset)}_history_comparison"
    fig.savefig(stem.with_suffix(".png"), dpi=260)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_isolated_cycles(dataset: str) -> None:
    cases = cases_for(dataset)
    fig, axes = plt.subplots(len(cases), 2, figsize=(11.5, 12.2), squeeze=False)
    for row, case in enumerate(cases):
        experiment, numerical = histories(case)
        for history, label, color, style, width in (
            (experiment, "Experiment", EXPERIMENT_COLOR, "-", 2.0),
            (numerical, "Corrected research successor", NUMERICAL_COLOR, "--", 1.8),
        ):
            mask = cycle_mask(history, case.comparison_cycle)
            if np.count_nonzero(mask) < 3:
                if label != "Experiment":
                    for axis in axes[row]:
                        axis.text(
                            0.98, 0.07,
                            f"Numerical cycle {case.comparison_cycle} unavailable\n"
                            f"(stopped at N={history['cycle'][-1]:.3g})",
                            color=NUMERICAL_COLOR, ha="right", va="bottom",
                            transform=axis.transAxes, fontsize=9,
                        )
                continue
            axes[row, 0].plot(
                history["gamma_percent"][mask], history["tau"][mask],
                color=color, ls=style, lw=width, label=label,
            )
            axes[row, 1].plot(
                history["vertical_effective_stress"][mask], history["tau"][mask],
                color=color, ls=style, lw=width,
            )
        axes[row, 0].set_ylabel(
            f"{case.case_id}, cycle {case.comparison_cycle}\n"
            r"$\tau$ (kPa)"
        )
        for axis in axes[row]:
            apply_grid(axis)
    axes[0, 0].set_title("Isolated stress–strain loop")
    axes[0, 1].set_title("Isolated effective-stress path")
    axes[-1, 0].set_xlabel(r"Shear strain, $\gamma$ (%)")
    axes[-1, 1].set_xlabel(r"Vertical effective stress, $\sigma'_v$ (kPa)")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{dataset}: cycle-matched loop and stress-path comparison")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.975))
    stem = OUTPUT / f"{dataset_slug(dataset)}_isolated_cycle_comparison"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_criterion_histories(dataset: str) -> None:
    cases = cases_for(dataset)
    metrics = read_metrics()
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6), squeeze=False)
    for axis, case in zip(axes.flat, cases):
        experiment, numerical = histories(case)
        for history, label, color, style, width in (
            (experiment, "Experiment", EXPERIMENT_COLOR, "-", 2.0),
            (numerical, "Corrected research successor", NUMERICAL_COLOR, "--", 1.8),
        ):
            cycle, measure = criterion_history(history, case.criterion)
            axis.plot(cycle, measure, color=color, ls=style, lw=width, label=label)
        axis.axhline(7.5, color=THRESHOLD_COLOR, lw=1.0, alpha=0.75)
        n_model = float(metrics[case.case_id]["N_model"])
        if np.isfinite(n_model):
            axis.scatter(
                [case.experimental_cycles, n_model], [7.5, 7.5],
                color=[EXPERIMENT_COLOR, NUMERICAL_COLOR], s=34, zorder=5,
            )
            n_text = f"Exp {case.experimental_cycles:.3g}; Num {n_model:.3g}"
        else:
            n_text = (
                f"Exp {case.experimental_cycles:.3g}; Num >"
                f"{float(metrics[case.case_id]['completed']):.3g} (runout)"
            )
        criterion_label = "DA" if case.criterion == "double_amplitude" else "SA"
        axis.set(
            title=f"{case.case_id}: {criterion_label} criterion\n{n_text}",
            xlabel="Cycle, N",
            ylabel=f"{criterion_label} shear strain (%)",
            ylim=(0.0, max(8.2, axis.get_ylim()[1])),
        )
        apply_grid(axis)
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{dataset}: strain accumulation and 7.5% criterion")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    stem = OUTPUT / f"{dataset_slug(dataset)}_strain_criterion_comparison"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_csr_n(dataset: str) -> None:
    cases = cases_for(dataset)
    metrics = read_metrics()
    fig, axis = plt.subplots(figsize=(8.6, 5.5))
    offsets = np.linspace(-0.0045, 0.0045, len(cases))
    numerical_labeled = False
    runout_labeled = False
    plotted_values = []
    for index, (case, offset) in enumerate(zip(cases, offsets)):
        n_exp = case.experimental_cycles
        n_model = float(metrics[case.case_id]["N_model"])
        plotted_values.append(n_exp)
        axis.scatter(
            case.csr + offset - 0.0015, n_exp,
            color=EXPERIMENT_COLOR, marker="o", s=58,
            label="Experiment" if index == 0 else None, zorder=4,
        )
        if np.isfinite(n_model):
            plotted_values.append(n_model)
            axis.scatter(
                case.csr + offset + 0.0015, n_model,
                color=NUMERICAL_COLOR, marker="s", s=54,
                label="Corrected research successor" if not numerical_labeled else None,
                zorder=4,
            )
            numerical_labeled = True
            axis.plot(
                [case.csr + offset - 0.0015, case.csr + offset + 0.0015],
                [n_exp, n_model], color="0.72", lw=0.9, zorder=1,
            )
            label_value = max(n_exp, n_model)
        else:
            completed = float(metrics[case.case_id]["completed"])
            plotted_values.append(completed)
            axis.scatter(
                case.csr + offset + 0.0015, completed,
                facecolors="none", edgecolors=NUMERICAL_COLOR,
                marker="^", s=72, linewidths=1.6,
                label="Numerical runout" if not runout_labeled else None,
                zorder=4,
            )
            runout_labeled = True
            axis.annotate(
                "", xy=(case.csr + offset + 0.0015, completed * 1.42),
                xytext=(case.csr + offset + 0.0015, completed),
                arrowprops=dict(arrowstyle="->", color=NUMERICAL_COLOR, lw=1.3),
            )
            label_value = max(n_exp, completed * 1.42)
        axis.annotate(
            case.case_id,
            xy=(case.csr + offset, label_value),
            xytext=(0, 9), textcoords="offset points",
            ha="center", va="bottom", fontsize=8.5,
        )
    axis.set_yscale("log")
    csr_values = [case.csr for case in cases]
    axis.set(
        xlabel="Cyclic stress ratio, CSR",
        ylabel="Cycles to criterion, N",
        xlim=(min(csr_values) - 0.035, max(csr_values) + 0.035),
        ylim=(min(plotted_values) / 1.7, max(plotted_values) * 1.85),
    )
    axis.set_title(
        f"{dataset}: cycles to the specified 7.5% strain criterion\n"
        "Individual tests; density and static bias are not constant",
        pad=12,
    )
    apply_grid(axis)
    axis.legend(frameon=False, fontsize=9, loc="best")
    fig.tight_layout()
    stem = OUTPUT / f"{dataset_slug(dataset)}_csr_n_comparison"
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def write_plot_index() -> None:
    metrics = read_metrics()
    with (OUTPUT / "README.md").open("w", encoding="utf-8") as stream:
        stream.write(
            "# Latest two-dataset comparison plots\n\n"
            "Numerical series use the private loose/unbiased correction successor "
            "from `loose_biased_shear_flow_final_full`. Experimental "
            "series are origin-shifted DSS records. Orange is experimental; blue "
            "is numerical.\n\n"
        )
        for dataset in DATASETS:
            stream.write(f"## {dataset}\n\n")
            for suffix in (
                "history_comparison", "isolated_cycle_comparison",
                "strain_criterion_comparison", "csr_n_comparison",
            ):
                stem = f"{dataset_slug(dataset)}_{suffix}"
                stream.write(f"- `{stem}.png` and `{stem}.pdf`\n")
            stream.write("\n")
        stream.write("## Cycle summary\n\n")
        stream.write("| Case | Dataset | CSR | Bias | Experiment N | Numerical N | Status |\n")
        stream.write("|---|---|---:|---:|---:|---:|---|\n")
        for case in audit.audit_module.CASES:
            record = metrics[case.case_id]
            value = float(record["N_model"])
            model_text = f"{value:.6g}" if np.isfinite(value) else "runout"
            stream.write(
                f"| {case.case_id} | {case.dataset} | {case.csr:g} | "
                f"{case.bias:g} | {case.experimental_cycles:.6g} | "
                f"{model_text} | {record['status']} |\n"
            )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for dataset in DATASETS:
        plot_history_matrix(dataset)
        plot_isolated_cycles(dataset)
        plot_criterion_histories(dataset)
        plot_csr_n(dataset)
    write_plot_index()
    print(OUTPUT)


if __name__ == "__main__":
    main()
