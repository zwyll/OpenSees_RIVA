"""State, restart, and integration checks for the biased-loop prototype."""

from __future__ import annotations

from dataclasses import asdict, fields
from pathlib import Path
import sys

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

from RIVASandBiasedLoop import (  # noqa: E402
    RIVASandBiasedLoopModel,
    RIVASandBiasedLoopParameters,
)
from rivasand_port.model import RIVASandModel, RIVASandState  # noqa: E402
from rivasand_port.reference import reference_parameters  # noqa: E402


def strain_increment(delta_gamma: float) -> np.ndarray:
    result = np.zeros((3, 3), dtype=float)
    result[0, 2] = result[2, 0] = 0.5 * delta_gamma
    return result


def assert_same_state(first: RIVASandState, second: RIVASandState) -> None:
    for item in fields(RIVASandState):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"state tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"state field {item.name} differs")


def initialized(model: RIVASandModel) -> RIVASandState:
    stress = np.diag([-19.4, -19.4, -40.0]).astype(float)
    stress[0, 2] = stress[2, 0] = 10.0
    model.initialize(stress, 0.536)
    reference = np.diag([-19.4, -19.4, -40.0]).astype(float)
    model.begin_cyclic_phase(reference_stress=reference)
    assert model.state is not None
    return model.state.copy()


def disabled_equivalence() -> None:
    base_parameters = reference_parameters()
    research_parameters = RIVASandBiasedLoopParameters(
        **asdict(base_parameters), branch_compliance_enabled=False
    )
    base = RIVASandModel(base_parameters)
    research = RIVASandBiasedLoopModel(research_parameters)
    base_state = initialized(base)
    research_state = initialized(research)
    assert_same_state(base_state, research_state)
    gamma = 0.0
    for step in range(1, 129):
        target = 0.005 * np.sin(2.0 * np.pi * step / 32.0)
        increment = strain_increment(target - gamma)
        gamma = target
        base_state, _ = base.advance_fixed(base_state, increment, 4)
        research_state, _ = research.advance_fixed(research_state, increment, 4)
        assert_same_state(base_state, research_state)


def run_enabled(nsub: int) -> tuple[RIVASandState, np.ndarray]:
    model = RIVASandBiasedLoopModel()
    state = initialized(model)
    gamma = 0.0
    stresses = []
    for step in range(1, 257):
        target = 0.0125 * np.sin(2.0 * np.pi * step / 32.0)
        increment = strain_increment(target - gamma)
        gamma = target
        state, _ = model.advance_fixed(state, increment, nsub)
        if not np.all(np.isfinite(state.stress)):
            raise AssertionError(f"nonfinite stress at step {step}")
        stresses.append(state.stress.copy())
    return state, np.asarray(stresses)


def restart_equivalence() -> None:
    model = RIVASandBiasedLoopModel()
    state = initialized(model)
    gamma = 0.0
    increments = []
    for step in range(1, 129):
        target = 0.01 * np.sin(2.0 * np.pi * step / 32.0)
        increments.append(strain_increment(target - gamma))
        gamma = target
    for increment in increments[:64]:
        state, _ = model.advance_fixed(state, increment, 4)
    restarted = state.copy()
    uninterrupted = state.copy()
    for increment in increments[64:]:
        restarted, _ = model.advance_fixed(restarted, increment, 4)
        uninterrupted, _ = model.advance_fixed(uninterrupted, increment, 4)
    assert_same_state(restarted, uninterrupted)


def substep_convergence() -> None:
    _, stress_two = run_enabled(2)
    _, stress_four = run_enabled(4)
    scale = max(float(np.max(np.abs(stress_four))), 1.0)
    relative = float(np.max(np.abs(stress_two - stress_four)) / scale)
    if relative > 0.08:
        raise AssertionError(f"2-versus-4 substep stress difference is {relative:.3g}")
    print(f"PASS: 2-versus-4 substep maximum stress difference = {relative:.3%}")


def main() -> None:
    disabled_equivalence()
    restart_equivalence()
    substep_convergence()
    print("PASS: disabled oracle equivalence, finite response, and restart")


if __name__ == "__main__":
    main()
