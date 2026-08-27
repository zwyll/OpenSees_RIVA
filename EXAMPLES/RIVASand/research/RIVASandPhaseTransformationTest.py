"""State, restart, disabled-equivalence, and integration checks for PT research."""

from __future__ import annotations

from dataclasses import fields, replace
from pathlib import Path
import sys

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

from RIVASandPhaseTransformation import (  # noqa: E402
    RIVASandPhaseTransformationModel,
    RIVASandPhaseTransformationParameters,
    RIVASandPhaseTransformationState,
)
from rivasand_port.model import (  # noqa: E402
    RIVASandModel,
    RIVASandParameters,
    RIVASandState,
)


def strain_increment(delta_gamma: float) -> np.ndarray:
    result = np.zeros((3, 3), dtype=float)
    result[0, 2] = result[2, 0] = 0.5 * delta_gamma
    return result


def assert_base_state_equal(first: RIVASandState, second: RIVASandState) -> None:
    for item in fields(RIVASandState):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"state tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"state field {item.name} differs")


def assert_pt_state_equal(
    first: RIVASandPhaseTransformationState,
    second: RIVASandPhaseTransformationState,
) -> None:
    for item in fields(RIVASandPhaseTransformationState):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"state tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"state field {item.name} differs")


def initialize(model):
    stress = np.diag([-19.4, -19.4, -40.0]).astype(float)
    stress[0, 2] = stress[2, 0] = 10.0
    model.initialize(stress, 0.536)
    reference = np.diag([-19.4, -19.4, -40.0]).astype(float)
    model.begin_cyclic_phase(reference_stress=reference)
    assert model.state is not None
    return model.state.copy()


def disabled_equivalence() -> None:
    parameters = RIVASandPhaseTransformationParameters()
    production_parameters = RIVASandParameters(**{
        item.name: getattr(parameters, item.name)
        for item in fields(RIVASandParameters)
    })
    phase_parameters = replace(parameters, phase_transformation_enabled=False)
    production = RIVASandModel(production_parameters)
    phase = RIVASandPhaseTransformationModel(phase_parameters)
    first = initialize(production)
    second = initialize(phase)
    assert_base_state_equal(first, second)
    gamma = 0.0
    for step in range(1, 129):
        target = 0.01 * np.sin(2.0 * np.pi * step / 32.0)
        increment = strain_increment(target - gamma)
        gamma = target
        first, _ = production.advance_fixed(first, increment, 4)
        second, _ = phase.advance_fixed(second, increment, 4)
        assert_base_state_equal(first, second)


def response(nsub: int) -> tuple[RIVASandPhaseTransformationState, np.ndarray]:
    model = RIVASandPhaseTransformationModel()
    state = initialize(model)
    gamma = 0.0
    stresses = []
    for step in range(1, 257):
        target = 0.0125 * np.sin(2.0 * np.pi * step / 32.0)
        increment = strain_increment(target - gamma)
        gamma = target
        state, _ = model.advance_fixed(state, increment, nsub)
        if not np.all(np.isfinite(state.stress)):
            raise AssertionError(f"nonfinite stress at step {step}")
        if abs(state.compatibility_residual) > 1.0e-12:
            raise AssertionError("volumetric compatibility residual is nonzero")
        stresses.append(state.stress.copy())
    return state, np.asarray(stresses)


def restart_equivalence() -> None:
    model = RIVASandPhaseTransformationModel()
    state = initialize(model)
    gamma = 0.0
    increments = []
    for step in range(1, 129):
        target = 0.01 * np.sin(2.0 * np.pi * step / 32.0)
        increments.append(strain_increment(target - gamma))
        gamma = target
    for increment in increments[:64]:
        state, _ = model.advance_fixed(state, increment, 4)
    first = state.copy()
    second = state.copy()
    for increment in increments[64:]:
        first, _ = model.advance_fixed(first, increment, 4)
        second, _ = model.advance_fixed(second, increment, 4)
    assert_pt_state_equal(first, second)


def substep_convergence() -> None:
    _, stress_two = response(2)
    _, stress_four = response(4)
    scale = max(float(np.max(np.abs(stress_four))), 1.0)
    relative = float(np.max(np.abs(stress_two - stress_four)) / scale)
    if relative > 0.10:
        raise AssertionError(f"2-versus-4 substep stress difference is {relative:.3g}")
    print(f"PASS: PT 2-versus-4 substep stress difference = {relative:.3%}")


def main() -> None:
    disabled_equivalence()
    restart_equivalence()
    substep_convergence()
    print("PASS: PT disabled equivalence, compatibility, finite response, and restart")


if __name__ == "__main__":
    main()
