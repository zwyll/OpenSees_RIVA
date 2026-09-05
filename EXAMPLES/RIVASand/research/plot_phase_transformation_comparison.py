"""Plot the restrained PT checkpoint against production and DSS experiments."""

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
RESEARCH = RESULTS / "phase_transformation_density_extended_final_full"
REDUCED_CONTRACTION = RESULTS / "pt_density12_cm050"
OUTPUT = RESULTS / "phase_transformation_density_extended_comparison"

EXPERIMENT_COLOR = "#D55E00"
PRODUCTION_COLOR = "#666666"
RESEARCH_COLOR = "#0072B2"
REDUCED_CONTRACTION_COLOR = "#009E73"
CASE_IDS = ("3484_b025", "4666_dense_b025", "4666_dense_b0375")
EFFECTIVE_PATH_CASE_IDS = CASE_IDS


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


def running_double_amplitude(history: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Return the driver criterion history, including the zero-strain start."""
    cycle = np.concatenate(([0.0], history["cycle"]))
    strain = np.concatenate(([0.0], history["gamma_percent"]))
    measure = np.maximum.accumulate(strain) - np.minimum.accumulate(strain)
    return cycle, measure


def histories_for(case):
    return (
        (
            audit.audit_module.shifted_experiment(case),
            "Experiment", EXPERIMENT_COLOR, "-", 1.6,
        ),
        (
            read_history(BASELINE / f"{case.case_id}.csv"),
            "Production RIVA-Sand", PRODUCTION_COLOR, ":", 1.4,
        ),
        (
            read_history(RESEARCH / f"{case.case_id}.csv"),
            "PT research", RESEARCH_COLOR, "--", 1.4,
        ),
    )


def plot_isolated_loops() -> None:
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for axis, case_id in zip(axes, CASE_IDS):
        case = cases[case_id]
        for history, label, color, style, width in histories_for(case):
            mask = cycle_mask(history, case.comparison_cycle)
            axis.plot(
                history["gamma_percent"][mask], history["tau"][mask],
                color=color, ls=style, lw=width, label=label,
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
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Biased DSS loops: restrained phase-transformation checkpoint")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    fig.savefig(OUTPUT / "phase_transformation_isolated_loops.png", dpi=240)
    fig.savefig(OUTPUT / "phase_transformation_isolated_loops.pdf")
    plt.close(fig)


def plot_histories() -> None:
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 10.4))
    for row, case_id in enumerate(CASE_IDS):
        case = cases[case_id]
        experiment = audit.audit_module.shifted_experiment(case)
        end = min(12.0, float(experiment["cycle"][-1]))
        for history, label, color, style, width in histories_for(case):
            mask = history["cycle"] <= min(end, float(history["cycle"][-1]))
            ru_key = "ru" if label == "Experiment" else "ru_vertical"
            axes[row, 0].plot(
                history["gamma_percent"][mask], history["tau"][mask],
                color=color, ls=style, lw=width, label=label,
            )
            axes[row, 1].plot(
                history["vertical_effective_stress"][mask],
                history["tau"][mask], color=color, ls=style, lw=width,
            )
            axes[row, 2].plot(
                history["cycle"][mask], history[ru_key][mask],
                color=color, ls=style, lw=width,
            )
        axes[row, 0].set_ylabel(f"{case_id}\nShear stress, τ (kPa)")
        for axis in axes[row]:
            axis.grid(color="0.88", lw=0.6)
    axes[0, 0].set_title("Stress–strain history")
    axes[0, 1].set_title("Effective-stress path")
    axes[0, 2].set_title(r"Pore-pressure ratio, $r_u$")
    axes[-1, 0].set_xlabel(r"Shear strain, $\gamma$ (%)")
    axes[-1, 1].set_xlabel(r"Vertical effective stress, $\sigma'_v$ (kPa)")
    axes[-1, 2].set_xlabel("Cycle, N")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle("Biased DSS histories: production and restrained PT research")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.975))
    fig.savefig(OUTPUT / "phase_transformation_histories.png", dpi=240)
    fig.savefig(OUTPUT / "phase_transformation_histories.pdf")
    plt.close(fig)


def plot_effective_paths() -> None:
    """Compare one complete experimental and numerical path loop."""
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.5))
    for axis, case_id in zip(axes, EFFECTIVE_PATH_CASE_IDS):
        case = cases[case_id]
        for history, label, color, style, width in histories_for(case):
            mask = cycle_mask(history, case.comparison_cycle)
            axis.plot(
                history["vertical_effective_stress"][mask],
                history["tau"][mask],
                color=color,
                ls=style,
                lw=width,
                label=label,
            )
        axis.set(
            title=(
                rf"$\alpha={case.bias:g}$, cycle {case.comparison_cycle}"
            ),
            xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
            ylabel=r"Shear stress, $\tau$ (kPa)",
        )
        axis.grid(color="0.88", lw=0.6)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Biased DSS: isolated effective-stress paths")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(OUTPUT / "phase_transformation_effective_paths.png", dpi=240)
    fig.savefig(OUTPUT / "phase_transformation_effective_paths.pdf")
    plt.close(fig)


def plot_3484_b025_effective_path() -> None:
    """Publication-scale view of the corrected intermediate-density path."""
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    case = cases["3484_b025"]
    series = histories_for(case)
    # Put the corrected result next to the experiment in the legend; retain
    # production as a subdued reference showing the original hourglass.
    order = (series[0], series[2], series[1])
    fig, axis = plt.subplots(figsize=(7.0, 5.8))
    for history, label, color, style, width in order:
        mask = cycle_mask(history, case.comparison_cycle)
        sigma_v = history["vertical_effective_stress"][mask]
        tau = history["tau"][mask]
        axis.plot(
            sigma_v, tau, color=color, ls=style,
            lw=2.3 if label != "Production RIVA-Sand" else 1.7,
            label=label,
        )
        axis.scatter(
            sigma_v[0], tau[0], s=25, facecolor=color, edgecolor="white",
            linewidth=0.6, zorder=4,
        )
    axis.set(
        title=(
            "PRJ-3484 b025: cycle-11 effective-stress path\n"
            r"$D_r=0.663$, $\alpha=0.25$, CSR$=0.20$"
        ),
        xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
        ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axis.grid(color="0.88", lw=0.6)
    axis.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT / "3484_b025_effective_stress_path.png", dpi=300)
    fig.savefig(OUTPUT / "3484_b025_effective_stress_path.pdf")
    plt.close(fig)


def plot_strain_criterion_explainer() -> None:
    """Show early triggering and why global contraction scaling is rejected."""
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    case25 = cases["3484_b025"]
    case30 = cases["3484_b030"]
    exp25 = audit.audit_module.shifted_experiment(case25)
    prod25 = read_history(BASELINE / "3484_b025.csv")
    pt25 = read_history(RESEARCH / "3484_b025.csv")
    reduced25 = read_history(REDUCED_CONTRACTION / "3484_b025.csv")
    exp30 = audit.audit_module.shifted_experiment(case30)
    prod30 = read_history(BASELINE / "3484_b030.csv")
    pt30 = read_history(RESEARCH / "3484_b030.csv")

    fig, axes = plt.subplots(1, 3, figsize=(15.4, 4.8))
    for history, label, color, style, width in (
        (exp25, "Experiment", EXPERIMENT_COLOR, "-", 2.0),
        (prod25, "Production", PRODUCTION_COLOR, ":", 1.6),
        (pt25, "PT research", RESEARCH_COLOR, "--", 2.0),
        (reduced25, "50% contraction", REDUCED_CONTRACTION_COLOR, "-.", 1.8),
    ):
        cycle, da = running_double_amplitude(history)
        mask = cycle <= 70.0
        axes[0].plot(cycle[mask], da[mask], color=color, ls=style, lw=width, label=label)
    axes[0].axhline(7.5, color="black", lw=1.0, alpha=0.7)
    axes[0].scatter(
        [66.2266, 27.15625, 26.21875, 48.1875], [7.5] * 4,
        color=[EXPERIMENT_COLOR, PRODUCTION_COLOR, RESEARCH_COLOR, REDUCED_CONTRACTION_COLOR],
        s=28, zorder=4,
    )
    axes[0].set(
        title="3484_b025: criterion timing",
        xlabel="Cycle, N", ylabel="Accumulated double-amplitude strain (%)",
        xlim=(0.0, 70.0), ylim=(0.0, 8.1),
    )
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")

    for history, label, color, style, width in (
        (exp25, "Experiment", EXPERIMENT_COLOR, "-", 2.0),
        (pt25, "PT research", RESEARCH_COLOR, "--", 2.0),
        (reduced25, "50% contraction", REDUCED_CONTRACTION_COLOR, "-.", 1.8),
    ):
        mask = cycle_mask(history, case25.comparison_cycle)
        axes[1].plot(
            history["gamma_percent"][mask], history["tau"][mask],
            color=color, ls=style, lw=width, label=label,
        )
    axes[1].set(
        title="3484_b025: cycle-11 loop",
        xlabel=r"Shear strain, $\gamma$ (%)", ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[1].legend(frameon=False, fontsize=8, loc="lower right")

    for history, label, color, style, width in (
        (exp30, "Experiment", EXPERIMENT_COLOR, "-", 2.0),
        (prod30, "Production", PRODUCTION_COLOR, ":", 1.6),
        (pt30, "PT research", RESEARCH_COLOR, "--", 2.0),
    ):
        cycle, da = running_double_amplitude(history)
        axes[2].plot(cycle, da, color=color, ls=style, lw=width, label=label)
    axes[2].axhline(7.5, color="black", lw=1.0, alpha=0.7)
    axes[2].scatter(
        [3.219, 1.0625, 1.0625], [7.5] * 3,
        color=[EXPERIMENT_COLOR, PRODUCTION_COLOR, RESEARCH_COLOR],
        s=28, zorder=4,
    )
    axes[2].set(
        title="3484_b030: criterion timing",
        xlabel="Cycle, N", ylabel="Accumulated double-amplitude strain (%)",
        xlim=(0.0, 3.5), ylim=(0.0, 8.1),
    )
    axes[2].legend(frameon=False, fontsize=8, loc="upper left")

    for axis in axes:
        axis.grid(color="0.88", lw=0.6)
    fig.suptitle("What premature strain triggering means")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(OUTPUT / "strain_criterion_explainer.png", dpi=300)
    fig.savefig(OUTPUT / "strain_criterion_explainer.pdf")
    plt.close(fig)


def write_effective_path_metrics() -> None:
    """Write the four cycle-path targets used in the density calibration."""
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    rows = []
    for case_id in EFFECTIVE_PATH_CASE_IDS:
        case = cases[case_id]
        for history, label, *_ in histories_for(case):
            mask = cycle_mask(history, case.comparison_cycle)
            sigma_v = history["vertical_effective_stress"][mask]
            tau = history["tau"][mask]
            area = 0.5 * abs(float(np.sum(
                sigma_v * np.roll(tau, -1)
                - np.roll(sigma_v, -1) * tau
            )))
            rows.append({
                "case_id": case_id,
                "series": label,
                "comparison_cycle": case.comparison_cycle,
                "mean_vertical_effective_stress": float(np.mean(sigma_v)),
                "vertical_effective_stress_range": float(np.ptp(sigma_v)),
                "enclosed_path_area": area,
                "tau_sigma_v_correlation": float(np.corrcoef(tau, sigma_v)[0, 1]),
            })
    with (OUTPUT / "effective_path_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_cycles() -> None:
    baseline = read_metrics(BASELINE / "metrics.csv")
    research = read_metrics(RESEARCH / "metrics.csv")
    cases = {case.case_id: case for case in audit.audit_module.CASES}
    x = np.arange(len(CASE_IDS), dtype=float)
    width = 0.24
    experiment = [cases[case_id].experimental_cycles for case_id in CASE_IDS]
    production = [float(baseline[case_id]["N_model"]) for case_id in CASE_IDS]
    candidate = [float(research[case_id]["N_model"]) for case_id in CASE_IDS]
    fig, axis = plt.subplots(figsize=(8.3, 4.5))
    axis.bar(x - width, experiment, width, color=EXPERIMENT_COLOR, label="Experiment")
    axis.bar(x, production, width, color=PRODUCTION_COLOR, label="Production")
    axis.bar(x + width, candidate, width, color=RESEARCH_COLOR, label="PT research")
    axis.set(
        ylabel="Cycles to strain criterion, N",
        xticks=x,
        xticklabels=CASE_IDS,
        title="Biased DSS cycles to strain criterion",
    )
    axis.grid(axis="y", color="0.88", lw=0.6)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT / "phase_transformation_cycles.png", dpi=240)
    fig.savefig(OUTPUT / "phase_transformation_cycles.pdf")
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_isolated_loops()
    plot_histories()
    plot_effective_paths()
    plot_3484_b025_effective_path()
    plot_strain_criterion_explainer()
    write_effective_path_metrics()
    plot_cycles()
    print(OUTPUT)


if __name__ == "__main__":
    main()
