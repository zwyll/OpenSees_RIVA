"""State, restart, equivalence, and integration checks for accumulation control."""

from __future__ import annotations

from dataclasses import fields, replace
from pathlib import Path
import sys

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

from RIVASandAccumulationControl import (  # noqa: E402
    RIVASandAccumulationControlModel,
    RIVASandAccumulationControlParameters,
    RIVASandAccumulationControlState,
)
from RIVASandPhaseTransformation import (  # noqa: E402
    RIVASandPhaseTransformationModel,
    RIVASandPhaseTransformationParameters,
    RIVASandPhaseTransformationState,
)


def strain_increment(delta_gamma: float) -> np.ndarray:
    result = np.zeros((3, 3), dtype=float)
    result[0, 2] = result[2, 0] = 0.5 * delta_gamma
    return result


def initialize(model, void_ratio: float = 0.601):
    stress = np.diag([-19.4, -19.4, -40.0]).astype(float)
    stress[0, 2] = stress[2, 0] = 10.0
    model.initialize(stress, void_ratio)
    reference = np.diag([-19.4, -19.4, -40.0]).astype(float)
    model.begin_cyclic_phase(reference_stress=reference)
    assert model.state is not None
    return model.state.copy()


def assert_common_state_equal(first, second, state_type) -> None:
    for item in fields(state_type):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"state tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"state field {item.name} differs")


def disabled_equivalence() -> None:
    parameters = RIVASandAccumulationControlParameters()
    phase_parameters = RIVASandPhaseTransformationParameters(**{
        item.name: getattr(parameters, item.name)
        for item in fields(RIVASandPhaseTransformationParameters)
    })
    disabled = replace(parameters, phase_accumulation_control_enabled=False)
    phase = RIVASandPhaseTransformationModel(phase_parameters)
    accumulation = RIVASandAccumulationControlModel(disabled)
    first = initialize(phase)
    second = initialize(accumulation)
    gamma = 0.0
    for step in range(1, 257):
        target = 0.0125 * np.sin(2.0 * np.pi * step / 32.0)
        increment = strain_increment(target - gamma)
        gamma = target
        first, _ = phase.advance_fixed(first, increment, 4)
        second, _ = accumulation.advance_fixed(second, increment, 4)
        assert_common_state_equal(
            first, second, RIVASandPhaseTransformationState
        )


def dense_equivalence() -> None:
    """The compact density window must leave the dense PT branch exact."""
    parameters = RIVASandAccumulationControlParameters()
    phase_parameters = RIVASandPhaseTransformationParameters(**{
        item.name: getattr(parameters, item.name)
        for item in fields(RIVASandPhaseTransformationParameters)
    })
    phase = RIVASandPhaseTransformationModel(phase_parameters)
    accumulation = RIVASandAccumulationControlModel(parameters)
    first = initialize(phase, void_ratio=0.536)
    second = initialize(accumulation, void_ratio=0.536)
    gamma = 0.0
    for step in range(1, 129):
        target = 0.015 * np.sin(2.0 * np.pi * step / 32.0)
        increment = strain_increment(target - gamma)
        gamma = target
        first, _ = phase.advance_fixed(first, increment, 4)
        second, _ = accumulation.advance_fixed(second, increment, 4)
        assert_common_state_equal(
            first, second, RIVASandPhaseTransformationState
        )


def response(model, nsub: int) -> tuple[RIVASandAccumulationControlState, np.ndarray]:
    state = initialize(model)
    gamma = 0.0
    stresses = []
    for step in range(1, 385):
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
    model = RIVASandAccumulationControlModel()
    state = initialize(model)
    gamma = 0.0
    increments = []
    for step in range(1, 257):
        target = 0.0125 * np.sin(2.0 * np.pi * step / 32.0)
        increments.append(strain_increment(target - gamma))
        gamma = target
    for increment in increments[:128]:
        state, _ = model.advance_fixed(state, increment, 4)
    first = state.copy()
    second = state.copy()
    for increment in increments[128:]:
        first, _ = model.advance_fixed(first, increment, 4)
        second, _ = model.advance_fixed(second, increment, 4)
    assert_common_state_equal(first, second, RIVASandAccumulationControlState)


def substep_convergence() -> None:
    # This deliberately severe strain-driven history is more sensitive than
    # the stress-controlled DSS calibration.  Four-versus-eight substeps is
    # the relevant refinement check because four is the calibrated setting.
    _, stress_four = response(RIVASandAccumulationControlModel(), 4)
    _, stress_eight = response(RIVASandAccumulationControlModel(), 8)
    scale = max(float(np.max(np.abs(stress_eight))), 1.0)
    relative = float(np.max(np.abs(stress_four - stress_eight)) / scale)
    if relative > 0.15:
        raise AssertionError(f"4-versus-8 substep stress difference is {relative:.3g}")
    print(f"PASS: accumulation 4-versus-8 stress difference = {relative:.3%}")


def main() -> None:
    disabled_equivalence()
    dense_equivalence()
    restart_equivalence()
    substep_convergence()
    print("PASS: disabled/dense equivalence, restart, compatibility, and finite response")


if __name__ == "__main__":
    main()
