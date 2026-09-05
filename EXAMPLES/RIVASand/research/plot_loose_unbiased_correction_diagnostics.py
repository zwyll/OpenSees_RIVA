"""Plot calibration diagnostics for the loose/unbiased research corrections."""

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

import RIVASandBaselineAudit as audit  # noqa: E402


RESULTS = HERE / "results"
OUTPUT = RESULTS / "loose_unbiased_correction_diagnostics"
EXPERIMENT_COLOR = "#D55E00"
PRIOR_COLOR = "#777777"
CORRECTED_COLOR = "#0072B2"
CANDIDATE_COLORS = ("#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#6A3D9A")


def read_history(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {
        name: np.atleast_1d(np.asarray(data[name], dtype=float))
        for name in data.dtype.names or ()
    }


def case(case_id: str):
    return next(item for item in audit.audit_module.CASES if item.case_id == case_id)


def cycle_mask(history: dict[str, np.ndarray], cycle: int) -> np.ndarray:
    return (history["cycle"] >= cycle) & (history["cycle"] < cycle + 1)


def apply_grid(axis) -> None:
    axis.grid(color="0.88", lw=0.6)
    axis.set_axisbelow(True)


def plot_loose_sensitivity() -> None:
    target = case("4666_loose_b015")
    experiment = audit.audit_module.shifted_experiment(target)
    exp_mask = cycle_mask(experiment, target.comparison_cycle)
    candidates = (
        (2.2, 4.0625),
        (2.4, 5.0625),
        (2.5, 7.09375),
        (2.6, 7.1875),
        (2.7, 7.0625),
    )
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5))
    axes[0].plot(
        experiment["gamma_percent"][exp_mask], experiment["tau"][exp_mask],
        color=EXPERIMENT_COLOR, lw=2.4, label="Experiment",
    )
    axes[1].plot(
        experiment["vertical_effective_stress"][exp_mask],
        experiment["tau"][exp_mask], color=EXPERIMENT_COLOR, lw=2.4,
    )
    axes[2].plot(
        experiment["cycle"], experiment["ru"],
        color=EXPERIMENT_COLOR, lw=2.4,
    )
    for (multiplier, n_model), color in zip(candidates, CANDIDATE_COLORS):
        history = read_history(
            RESULTS / f"loose_hard_{int(round(multiplier * 10)):02d}"
            / f"{target.case_id}.csv"
        )
        mask = cycle_mask(history, target.comparison_cycle)
        label = rf"$H$ multiplier {multiplier:g}; $N={n_model:g}$"
        axes[0].plot(
            history["gamma_percent"][mask], history["tau"][mask],
            color=color, ls="--", lw=1.5, label=label,
        )
        axes[1].plot(
            history["vertical_effective_stress"][mask], history["tau"][mask],
            color=color, ls="--", lw=1.5,
        )
        axes[2].plot(
            history["cycle"], history["ru_vertical"],
            color=color, ls="--", lw=1.5,
        )
    axes[0].set(
        title="Cycle-3 stress–strain loop",
        xlabel=r"Shear strain, $\gamma$ (%)", ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[1].set(
        title="Cycle-3 effective-stress path",
        xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
        ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[2].set(
        title=r"Pore-pressure ratio, $r_u$",
        xlabel="Cycle, N", ylabel=r"$r_u$", xlim=(0.0, 9.5),
    )
    for axis in axes:
        apply_grid(axis)
    axes[0].legend(frameon=False, fontsize=8.2)
    fig.suptitle("PRJ-4666 loose biased DSS: isolated hardening calibration")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    stem = OUTPUT / "loose_hardening_sensitivity"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_targeted_before_after() -> None:
    specifications = (
        (
            case("4666_loose_b015"),
            RESULTS / "loose_phase_repl_00"
            / "4666_loose_b015.csv",
            RESULTS / "loose_unbiased_correction_final_full"
            / "4666_loose_b015.csv",
        ),
        (
            case("4666_dense_u038"),
            RESULTS / "accumulation_control_host_committed_final_full"
            / "4666_dense_u038.csv",
            RESULTS / "unbiased_full_s07_c08" / "4666_dense_u038.csv",
        ),
    )
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.4), squeeze=False)
    for row, (target, prior_path, corrected_path) in enumerate(specifications):
        experiment = audit.audit_module.shifted_experiment(target)
        prior = read_history(prior_path)
        corrected = read_history(corrected_path)
        for history, ru_key, label, color, style, width in (
            (experiment, "ru", "Experiment", EXPERIMENT_COLOR, "-", 2.3),
            (prior, "ru_vertical", "Prior research checkpoint", PRIOR_COLOR, ":", 1.8),
            (corrected, "ru_vertical", "Corrected research checkpoint", CORRECTED_COLOR, "--", 2.0),
        ):
            mask = cycle_mask(history, target.comparison_cycle)
            if np.count_nonzero(mask) >= 3:
                axes[row, 0].plot(
                    history["gamma_percent"][mask], history["tau"][mask],
                    color=color, ls=style, lw=width, label=label,
                )
                axes[row, 1].plot(
                    history["vertical_effective_stress"][mask], history["tau"][mask],
                    color=color, ls=style, lw=width,
                )
            axes[row, 2].plot(
                history["cycle"], history[ru_key],
                color=color, ls=style, lw=width,
            )
        axes[row, 0].set_ylabel(
            f"{target.case_id}, cycle {target.comparison_cycle}\n"
            r"Shear stress, $\tau$ (kPa)"
        )
        for axis in axes[row]:
            apply_grid(axis)
    axes[0, 0].set_title("Isolated stress–strain loop")
    axes[0, 1].set_title("Isolated effective-stress path")
    axes[0, 2].set_title(r"Pore-pressure ratio, $r_u$")
    axes[-1, 0].set_xlabel(r"Shear strain, $\gamma$ (%)")
    axes[-1, 1].set_xlabel(r"Vertical effective stress, $\sigma'_v$ (kPa)")
    axes[-1, 2].set_xlabel("Cycle, N")
    axes[0, 0].legend(frameon=False, fontsize=8.6)
    fig.suptitle("Targeted loose-sand triggering and dense zero-bias path corrections")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    stem = OUTPUT / "targeted_before_after"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_loose_phase_replacement() -> None:
    target = case("4666_loose_b015")
    experiment = audit.audit_module.shifted_experiment(target)
    exp_mask = cycle_mask(experiment, target.comparison_cycle)
    candidates = (
        ("00", "0.00", "#777777"),
        ("025", "0.25", "#009E73"),
        ("050", "0.50", "#56B4E9"),
        ("075", "0.75", "#CC79A7"),
        ("10", "1.00", "#0072B2"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.5))
    for axis_index, (x_key, y_key) in enumerate((
        ("gamma_percent", "tau"),
        ("vertical_effective_stress", "tau"),
        ("cycle", "ru"),
    )):
        mask = exp_mask if axis_index < 2 else np.ones_like(experiment["cycle"], dtype=bool)
        axes[axis_index].plot(
            experiment[x_key][mask], experiment[y_key][mask],
            color=EXPERIMENT_COLOR, lw=2.4, label="Experiment",
        )
    for suffix, label, color in candidates:
        history = read_history(
            RESULTS / f"loose_phase_repl_{suffix}" / f"{target.case_id}.csv"
        )
        mask = cycle_mask(history, target.comparison_cycle)
        axes[0].plot(
            history["gamma_percent"][mask], history["tau"][mask],
            color=color, ls="--", lw=1.5,
            label=f"Legacy-wave replacement {label}",
        )
        axes[1].plot(
            history["vertical_effective_stress"][mask], history["tau"][mask],
            color=color, ls="--", lw=1.5,
        )
        axes[2].plot(
            history["cycle"], history["ru_vertical"],
            color=color, ls="--", lw=1.5,
        )
    axes[0].set(
        title="Cycle-3 stress–strain loop",
        xlabel=r"Shear strain, $\gamma$ (%)", ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[1].set(
        title="Cycle-3 effective-stress path",
        xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
        ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[2].set(
        title=r"Pore-pressure ratio, $r_u$",
        xlabel="Cycle, N", ylabel=r"$r_u$", xlim=(0.0, 8.0),
    )
    for axis in axes:
        apply_grid(axis)
    axes[0].legend(frameon=False, fontsize=8.0)
    fig.suptitle("PRJ-4666 loose biased DSS: legacy pressure-wave replacement")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    stem = OUTPUT / "loose_phase_replacement"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_loose_phase_activity() -> None:
    target = case("4666_loose_b015")
    experiment = audit.audit_module.shifted_experiment(target)
    exp_mask = cycle_mask(experiment, target.comparison_cycle)
    candidates = (
        ("002", "0.02", "#777777"),
        ("005", "0.05", "#009E73"),
        ("010", "0.10", "#56B4E9"),
        ("020", "0.20", "#0072B2"),
        ("035", "0.35", "#CC79A7"),
        ("050", "0.50", "#6A3D9A"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.5))
    axes[0].plot(
        experiment["gamma_percent"][exp_mask], experiment["tau"][exp_mask],
        color=EXPERIMENT_COLOR, lw=2.4, label="Experiment",
    )
    axes[1].plot(
        experiment["vertical_effective_stress"][exp_mask],
        experiment["tau"][exp_mask], color=EXPERIMENT_COLOR, lw=2.4,
    )
    axes[2].plot(
        experiment["cycle"], experiment["ru"],
        color=EXPERIMENT_COLOR, lw=2.4,
    )
    for suffix, label, color in candidates:
        history = read_history(
            RESULTS / f"loose_phase_act_{suffix}" / f"{target.case_id}.csv"
        )
        mask = cycle_mask(history, target.comparison_cycle)
        axes[0].plot(
            history["gamma_percent"][mask], history["tau"][mask],
            color=color, ls="--", lw=1.5, label=f"PT activity {label}",
        )
        axes[1].plot(
            history["vertical_effective_stress"][mask], history["tau"][mask],
            color=color, ls="--", lw=1.5,
        )
        axes[2].plot(
            history["cycle"], history["ru_vertical"],
            color=color, ls="--", lw=1.5,
        )
    axes[0].set(
        title="Cycle-3 stress–strain loop",
        xlabel=r"Shear strain, $\gamma$ (%)", ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[1].set(
        title="Cycle-3 effective-stress path",
        xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
        ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[2].set(
        title=r"Pore-pressure ratio, $r_u$",
        xlabel="Cycle, N", ylabel=r"$r_u$", xlim=(0.0, 8.0),
    )
    for axis in axes:
        apply_grid(axis)
    axes[0].legend(frameon=False, fontsize=8.0)
    fig.suptitle("PRJ-4666 loose biased DSS: continuous signed PT activity")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    stem = OUTPUT / "loose_phase_activity"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_loose_phase_finalists() -> None:
    target = case("4666_loose_b015")
    experiment = audit.audit_module.shifted_experiment(target)
    exp_mask = cycle_mask(experiment, target.comparison_cycle)
    candidates = (
        (
            RESULTS / "loose_phase_repl_00" / f"{target.case_id}.csv",
            "Prior corrected timing", "#777777", ":",
        ),
        (
            RESULTS / "loose_phase_grid_w030_m110" / f"{target.case_id}.csv",
            "PT wave 0.30, mean 1.10", "#0072B2", "--",
        ),
        (
            RESULTS / "loose_phase_grid_w040_m140" / f"{target.case_id}.csv",
            "PT wave 0.40, mean 1.40", "#009E73", "-.",
        ),
        (
            RESULTS / "loose_phase_grid_w030_m140" / f"{target.case_id}.csv",
            "PT wave 0.30, mean 1.40", "#CC79A7", "--",
        ),
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.5))
    axes[0].plot(
        experiment["gamma_percent"][exp_mask], experiment["tau"][exp_mask],
        color=EXPERIMENT_COLOR, lw=2.5, label="Experiment",
    )
    axes[1].plot(
        experiment["vertical_effective_stress"][exp_mask],
        experiment["tau"][exp_mask], color=EXPERIMENT_COLOR, lw=2.5,
    )
    axes[2].plot(
        experiment["cycle"], experiment["ru"],
        color=EXPERIMENT_COLOR, lw=2.5,
    )
    for path, label, color, style in candidates:
        history = read_history(path)
        mask = cycle_mask(history, target.comparison_cycle)
        axes[0].plot(
            history["gamma_percent"][mask], history["tau"][mask],
            color=color, ls=style, lw=1.8, label=label,
        )
        axes[1].plot(
            history["vertical_effective_stress"][mask], history["tau"][mask],
            color=color, ls=style, lw=1.8,
        )
        axes[2].plot(
            history["cycle"], history["ru_vertical"],
            color=color, ls=style, lw=1.8,
        )
    axes[0].set(
        title="Cycle-3 stress–strain loop",
        xlabel=r"Shear strain, $\gamma$ (%)", ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[1].set(
        title="Cycle-3 effective-stress path",
        xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
        ylabel=r"Shear stress, $\tau$ (kPa)",
    )
    axes[2].set(
        title=r"Pore-pressure ratio, $r_u$",
        xlabel="Cycle, N", ylabel=r"$r_u$", xlim=(0.0, 8.0),
    )
    for axis in axes:
        apply_grid(axis)
    axes[0].legend(frameon=False, fontsize=8.0)
    fig.suptitle("PRJ-4666 loose biased DSS: signed PT finalists")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    stem = OUTPUT / "loose_phase_finalists"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_loose_sensitivity()
    plot_targeted_before_after()
    plot_loose_phase_replacement()
    plot_loose_phase_activity()
    plot_loose_phase_finalists()
    print(OUTPUT)


if __name__ == "__main__":
    main()
