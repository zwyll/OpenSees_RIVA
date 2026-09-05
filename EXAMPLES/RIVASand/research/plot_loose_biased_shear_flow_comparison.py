"""Plot the loose-biased shear-flow checkpoint against its parent and DSS data."""

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
OUTPUT = RESULTS / "loose_biased_shear_flow_comparison"
EXPERIMENT_COLOR = "#D55E00"
PARENT_COLOR = "#777777"
SUCCESSOR_COLOR = "#0072B2"


def read_history(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {
        name: np.atleast_1d(np.asarray(data[name], dtype=float))
        for name in data.dtype.names or ()
    }


def cycle_mask(history: dict[str, np.ndarray], cycle: int) -> np.ndarray:
    return (history["cycle"] >= cycle) & (history["cycle"] < cycle + 1)


def apply_grid(axis) -> None:
    axis.grid(color="0.88", lw=0.6)
    axis.set_axisbelow(True)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    target = next(
        item for item in audit.audit_module.CASES
        if item.case_id == "4666_loose_b015"
    )
    experiment = audit.audit_module.shifted_experiment(target)
    parent = read_history(
        RESULTS / "loose_unbiased_correction_final_full"
        / "4666_loose_b015.csv"
    )
    successor = read_history(
        RESULTS / "loose_biased_shear_flow_final_full"
        / "4666_loose_b015.csv"
    )
    histories = (
        (experiment, "ru", "Experiment", EXPERIMENT_COLOR, "-", 2.4),
        (parent, "ru_vertical", "Signed-PT parent", PARENT_COLOR, ":", 1.9),
        (
            successor, "ru_vertical", "Loose shear-flow successor",
            SUCCESSOR_COLOR, "--", 2.1,
        ),
    )
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.6), squeeze=False)
    for history, ru_key, label, color, style, width in histories:
        mask = cycle_mask(history, target.comparison_cycle)
        axes[0, 0].plot(
            history["gamma_percent"][mask], history["tau"][mask],
            color=color, ls=style, lw=width, label=label,
        )
        axes[0, 1].plot(
            history["vertical_effective_stress"][mask], history["tau"][mask],
            color=color, ls=style, lw=width,
        )
        axes[0, 2].plot(
            history["cycle"], history[ru_key], color=color, ls=style, lw=width,
        )
        full = history["cycle"] <= 8.0
        axes[1, 0].plot(
            history["gamma_percent"][full], history["tau"][full],
            color=color, ls=style, lw=width,
        )
        axes[1, 1].plot(
            history["vertical_effective_stress"][full], history["tau"][full],
            color=color, ls=style, lw=width,
        )
        axes[1, 2].plot(
            history["cycle"][full], history[ru_key][full],
            color=color, ls=style, lw=width,
        )
    axes[0, 0].set_title("Cycle-3 stress–strain loop")
    axes[0, 1].set_title("Cycle-3 effective-stress path")
    axes[0, 2].set_title(r"Pore-pressure ratio, $r_u$")
    axes[1, 0].set_title("First eight stress–strain cycles")
    axes[1, 1].set_title("First eight effective-stress cycles")
    axes[1, 2].set_title(r"First eight $r_u$ cycles")
    for axis in axes[:, 0]:
        axis.set(xlabel=r"Shear strain, $\gamma$ (%)", ylabel=r"$\tau$ (kPa)")
    for axis in axes[:, 1]:
        axis.set(
            xlabel=r"Vertical effective stress, $\sigma'_v$ (kPa)",
            ylabel=r"$\tau$ (kPa)",
        )
    for axis in axes[:, 2]:
        axis.set(xlabel="Cycle, N", ylabel=r"$r_u$", xlim=(0.0, 8.0))
    for axis in axes.flat:
        apply_grid(axis)
    axes[0, 0].legend(frameon=False, fontsize=8.7)
    fig.suptitle("PRJ-4666 loose biased DSS: shear-flow successor")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    stem = OUTPUT / "loose_biased_shear_flow_comparison"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()
