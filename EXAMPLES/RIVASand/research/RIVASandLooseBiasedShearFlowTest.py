"""Regression checks for the loose-biased shear-flow research successor."""

from __future__ import annotations

from dataclasses import fields, replace
from pathlib import Path
import sys

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

from RIVASandLooseBiasedShearFlow import (  # noqa: E402
    RIVASandLooseBiasedShearFlowModel,
    RIVASandLooseBiasedShearFlowParameters,
    RIVASandLooseBiasedShearFlowState,
)
from RIVASandLooseUnbiasedCorrection import (  # noqa: E402
    RIVASandLooseUnbiasedCorrectionModel,
    RIVASandLooseUnbiasedCorrectionParameters,
    RIVASandLooseUnbiasedCorrectionState,
)


def strain_increment(delta_gamma: float) -> np.ndarray:
    result = np.zeros((3, 3), dtype=float)
    result[0, 2] = result[2, 0] = 0.5 * delta_gamma
    return result


def initialize(model, *, void_ratio: float, bias_shear: float):
    reference = np.diag([-19.4, -19.4, -40.0]).astype(float)
    stress = reference.copy()
    stress[0, 2] = stress[2, 0] = bias_shear
    model.initialize(stress, void_ratio)
    model.begin_cyclic_phase(reference_stress=reference)
    assert model.state is not None
    return model.state.copy()


def parent_parameters(
    parameters: RIVASandLooseBiasedShearFlowParameters,
) -> RIVASandLooseUnbiasedCorrectionParameters:
    return RIVASandLooseUnbiasedCorrectionParameters(**{
        item.name: getattr(parameters, item.name)
        for item in fields(RIVASandLooseUnbiasedCorrectionParameters)
    })


def assert_common_state_equal(first, second, state_type) -> None:
    for item in fields(state_type):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"state tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"state field {item.name} differs")


def increments(amplitude: float = 0.015) -> list[np.ndarray]:
    result = []
    gamma = 0.0
    for step in range(1, 161):
        target = amplitude * np.sin(2.0 * np.pi * step / 32.0)
        result.append(strain_increment(target - gamma))
        gamma = target
    return result


def disabled_and_out_of_window_equivalence() -> None:
    parameters = RIVASandLooseBiasedShearFlowParameters()
    histories = (
        (replace(parameters, loose_shear_flow_enabled=False), 0.657, 6.0),
        (parameters, 0.601, 10.0),
        (parameters, 0.549, 0.0),
    )
    for candidate_parameters, void_ratio, bias_shear in histories:
        parent = RIVASandLooseUnbiasedCorrectionModel(
            parent_parameters(candidate_parameters)
        )
        successor = RIVASandLooseBiasedShearFlowModel(candidate_parameters)
        left = initialize(parent, void_ratio=void_ratio, bias_shear=bias_shear)
        right = initialize(successor, void_ratio=void_ratio, bias_shear=bias_shear)
        for increment in increments()[:96]:
            left, _ = parent.advance_fixed(left, increment, 4)
            right, _ = successor.advance_fixed(right, increment, 4)
            assert_common_state_equal(
                left, right, RIVASandLooseUnbiasedCorrectionState
            )


def finite_restart_and_committed_memory() -> None:
    model = RIVASandLooseBiasedShearFlowModel()
    state = initialize(model, void_ratio=0.657, bias_shear=6.0)
    history = increments(0.030)
    for increment in history[:96]:
        state, _ = model.advance_fixed(state, increment, 4)
        if not np.all(np.isfinite(state.stress)):
            raise AssertionError("nonfinite loose shear-flow response")
        if abs(state.compatibility_residual) > 1.0e-12:
            raise AssertionError("loose shear-flow compatibility residual is nonzero")
        if not 0.0 <= state.loose_shear_hardening_state <= 1.0:
            raise AssertionError("loose shear hardening state left [0,1]")
    if state.loose_shear_hardening_state <= 0.0:
        raise AssertionError("loose shear hardening memory did not evolve")
    first = state.copy()
    second = state.copy()
    for increment in history[96:]:
        first, _ = model.advance_fixed(first, increment, 4)
        second, _ = model.advance_fixed(second, increment, 4)
    assert_common_state_equal(first, second, RIVASandLooseBiasedShearFlowState)


def substep_refinement() -> float:
    first_model = RIVASandLooseBiasedShearFlowModel()
    second_model = RIVASandLooseBiasedShearFlowModel()
    first = initialize(first_model, void_ratio=0.657, bias_shear=6.0)
    second = initialize(second_model, void_ratio=0.657, bias_shear=6.0)
    maximum = 0.0
    scale = 1.0
    for increment in increments()[:96]:
        first, _ = first_model.advance_fixed(first, increment, 2)
        second, _ = second_model.advance_fixed(second, increment, 4)
        maximum = max(maximum, float(np.linalg.norm(first.stress - second.stress)))
        scale = max(scale, float(np.linalg.norm(second.stress)))
    relative = maximum / scale
    if relative > 0.15:
        raise AssertionError(f"2-versus-4 stress difference is {relative:.3%}")
    return relative


def main() -> None:
    disabled_and_out_of_window_equivalence()
    finite_restart_and_committed_memory()
    relative = substep_refinement()
    print(
        "PASS: disabled/out-of-window equivalence, finite response, restart, "
        f"compatibility, committed memory, and 2-versus-4 stress={relative:.3%}"
    )


if __name__ == "__main__":
    main()
