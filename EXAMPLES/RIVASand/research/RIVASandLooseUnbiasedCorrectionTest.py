"""Regression checks for the loose/unbiased RIVA-Sand research successor."""

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


def assert_common_state_equal(first, second, state_type) -> None:
    for item in fields(state_type):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"state tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"state field {item.name} differs")


def parent_parameters(
    parameters: RIVASandLooseUnbiasedCorrectionParameters,
) -> RIVASandAccumulationControlParameters:
    return RIVASandAccumulationControlParameters(**{
        item.name: getattr(parameters, item.name)
        for item in fields(RIVASandAccumulationControlParameters)
    })


def run_pair(
    first,
    second,
    *,
    void_ratio: float,
    bias_shear: float,
    amplitude: float = 0.0125,
) -> None:
    left = initialize(first, void_ratio=void_ratio, bias_shear=bias_shear)
    right = initialize(second, void_ratio=void_ratio, bias_shear=bias_shear)
    gamma = 0.0
    for step in range(1, 129):
        target = amplitude * np.sin(2.0 * np.pi * step / 32.0)
        increment = strain_increment(target - gamma)
        gamma = target
        left, _ = first.advance_fixed(left, increment, 4)
        right, _ = second.advance_fixed(right, increment, 4)
        assert_common_state_equal(left, right, RIVASandAccumulationControlState)


def disabled_equivalence() -> None:
    parameters = RIVASandLooseUnbiasedCorrectionParameters()
    disabled = replace(
        parameters,
        loose_stabilization_enabled=False,
        unbiased_phase_enabled=False,
    )
    parent = RIVASandAccumulationControlModel(parent_parameters(parameters))
    successor = RIVASandLooseUnbiasedCorrectionModel(disabled)
    run_pair(parent, successor, void_ratio=0.601, bias_shear=10.0)


def out_of_window_equivalence() -> None:
    """Medium-density and dense-biased paths must inherit the parent exactly."""
    parameters = RIVASandLooseUnbiasedCorrectionParameters()
    for void_ratio, bias_shear in ((0.601, 10.0), (0.549, 10.0)):
        parent = RIVASandAccumulationControlModel(parent_parameters(parameters))
        successor = RIVASandLooseUnbiasedCorrectionModel(parameters)
        run_pair(
            parent,
            successor,
            void_ratio=void_ratio,
            bias_shear=bias_shear,
        )


def fixed_direction_and_restart() -> None:
    model = RIVASandLooseUnbiasedCorrectionModel()
    state = initialize(model, void_ratio=0.549, bias_shear=0.0)
    gamma = 0.0
    increments = []
    for step in range(1, 129):
        target = 0.015 * np.sin(2.0 * np.pi * step / 32.0)
        increments.append(strain_increment(target - gamma))
        gamma = target
    direction = None
    for increment in increments[:64]:
        state, _ = model.advance_fixed(state, increment, 4)
        current_norm = np.linalg.norm(state.unbiased_phase_direction)
        if direction is None and current_norm > 0.0:
            direction = state.unbiased_phase_direction.copy()
        if direction is not None:
            alignment = float(np.sum(state.unbiased_phase_direction * direction))
            if alignment < 1.0 - 1.0e-13:
                raise AssertionError("fixed zero-bias phase direction reversed")
        if not np.all(np.isfinite(state.stress)):
            raise AssertionError("nonfinite dense zero-bias response")
        if abs(state.compatibility_residual) > 1.0e-12:
            raise AssertionError("volumetric compatibility residual is nonzero")
    if direction is None or not np.isclose(np.linalg.norm(direction), 1.0):
        raise AssertionError("dense zero-bias phase direction was not captured")
    first = state.copy()
    second = state.copy()
    for increment in increments[64:]:
        first, _ = model.advance_fixed(first, increment, 4)
        second, _ = model.advance_fixed(second, increment, 4)
    assert_common_state_equal(first, second, RIVASandLooseUnbiasedCorrectionState)


def gate_separation() -> None:
    model = RIVASandLooseUnbiasedCorrectionModel()
    loose = initialize(model, void_ratio=0.657, bias_shear=6.0)
    if model.loose_stabilization_gate(loose) < 0.99:
        raise AssertionError("loose biased gate did not activate")
    if model.unbiased_phase_gate(loose) != 0.0:
        raise AssertionError("unbiased phase gate activated in loose biased sand")
    if model.loose_phase_gate(loose) < 0.99:
        raise AssertionError("loose signed phase-volume gate did not activate")
    dense = initialize(model, void_ratio=0.549, bias_shear=0.0)
    if model.unbiased_phase_gate(dense) <= 0.0:
        raise AssertionError("dense zero-bias gate did not activate")
    if model.loose_stabilization_gate(dense) != 0.0:
        raise AssertionError("loose gate activated in dense sand")
    if model.loose_phase_gate(dense) != 0.0:
        raise AssertionError("loose phase-volume gate activated in dense sand")


def loose_phase_restart() -> None:
    model = RIVASandLooseUnbiasedCorrectionModel()
    state = initialize(model, void_ratio=0.657, bias_shear=6.0)
    gamma = 0.0
    increments = []
    for step in range(1, 129):
        target = 0.015 * np.sin(2.0 * np.pi * step / 32.0)
        increments.append(strain_increment(target - gamma))
        gamma = target
    for increment in increments[:64]:
        state, _ = model.advance_fixed(state, increment, 4)
        if not np.all(np.isfinite(state.stress)):
            raise AssertionError("nonfinite loose phase-volume response")
        if abs(state.compatibility_residual) > 1.0e-12:
            raise AssertionError("loose phase-volume compatibility residual is nonzero")
    if abs(state.phase_reversible_volume) <= 1.0e-12:
        raise AssertionError("loose signed phase volume did not evolve")
    first = state.copy()
    second = state.copy()
    for increment in increments[64:]:
        first, _ = model.advance_fixed(first, increment, 4)
        second, _ = model.advance_fixed(second, increment, 4)
    assert_common_state_equal(first, second, RIVASandLooseUnbiasedCorrectionState)


def main() -> None:
    disabled_equivalence()
    out_of_window_equivalence()
    fixed_direction_and_restart()
    gate_separation()
    loose_phase_restart()
    print(
        "PASS: disabled/out-of-window equivalence, gate separation, fixed "
        "direction, loose/dense restart, compatibility, and finite response"
    )


if __name__ == "__main__":
    main()
