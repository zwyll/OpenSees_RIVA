"""Numerical gates for the directional mapping-backstress successor."""

from __future__ import annotations

from dataclasses import fields, replace
from pathlib import Path
import sys

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

from RIVASandMappingBackstress import (  # noqa: E402
    RIVASandMappingBackstressModel,
    RIVASandMappingBackstressParameters,
    RIVASandMappingBackstressState,
)
from RIVASandLooseBiasedShearFlow import (  # noqa: E402
    RIVASandLooseBiasedShearFlowModel,
    RIVASandLooseBiasedShearFlowParameters,
    RIVASandLooseBiasedShearFlowState,
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
    parameters: RIVASandMappingBackstressParameters,
) -> RIVASandLooseBiasedShearFlowParameters:
    return RIVASandLooseBiasedShearFlowParameters(**{
        item.name: getattr(parameters, item.name)
        for item in fields(RIVASandLooseBiasedShearFlowParameters)
    })


def assert_parent_state_equal(first, second) -> None:
    for item in fields(RIVASandLooseBiasedShearFlowState):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"state tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"state field {item.name} differs")


def sinusoidal_increments(amplitude: float, points: int = 32, cycles: int = 4):
    result = []
    gamma = 0.0
    for step in range(1, points * cycles + 1):
        target = amplitude * np.sin(2.0 * np.pi * step / points)
        result.append(strain_increment(target - gamma))
        gamma = target
    return result


def disabled_and_preserved_window_equivalence() -> None:
    defaults = RIVASandMappingBackstressParameters()
    histories = (
        (replace(defaults, mapping_backstress_enabled=False), 0.660, 10.0),
        (defaults, 0.643, 6.0),
        (defaults, 0.601, 10.0),
        (defaults, 0.536, 10.0),
    )
    for candidate, void_ratio, bias_shear in histories:
        parent = RIVASandLooseBiasedShearFlowModel(parent_parameters(candidate))
        successor = RIVASandMappingBackstressModel(candidate)
        left = initialize(parent, void_ratio=void_ratio, bias_shear=bias_shear)
        right = initialize(successor, void_ratio=void_ratio, bias_shear=bias_shear)
        for increment in sinusoidal_increments(0.025, cycles=3):
            left, _ = parent.advance_fixed(left, increment, 4)
            right, _ = successor.advance_fixed(right, increment, 4)
            assert_parent_state_equal(left, right)


def finite_restart_compatibility_and_outer_surface() -> None:
    model = RIVASandMappingBackstressModel()
    state = initialize(model, void_ratio=0.660, bias_shear=10.0)
    history = sinusoidal_increments(0.035, cycles=6)
    activated = False
    for increment in history[:96]:
        state, _ = model.advance_fixed(state, increment, 4)
        activated = activated or np.linalg.norm(state.mapping_backstress) > 0.0
        if not np.all(np.isfinite(state.stress)):
            raise AssertionError("nonfinite mapping-backstress response")
        if abs(state.compatibility_residual) > 1.0e-12:
            raise AssertionError("mapping-backstress compatibility residual is nonzero")
        if state.mapping_outer_residual > 2.0e-8:
            raise AssertionError(
                f"outer-surface residual={state.mapping_outer_residual:.3e}"
            )
    if not activated:
        raise AssertionError("mapping backstress did not evolve")
    first = state.copy()
    second = state.copy()
    for increment in history[96:]:
        first, _ = model.advance_fixed(first, increment, 4)
        second, _ = model.advance_fixed(second, increment, 4)
    for item in fields(RIVASandMappingBackstressState):
        left, right = getattr(first, item.name), getattr(second, item.name)
        if isinstance(left, np.ndarray):
            if not np.array_equal(left, right):
                raise AssertionError(f"restart tensor {item.name} differs")
        elif left != right:
            raise AssertionError(f"restart field {item.name} differs")


def monotone_local_stress_map() -> None:
    model = RIVASandMappingBackstressModel()
    state = initialize(model, void_ratio=0.660, bias_shear=10.0)
    for increment in sinusoidal_increments(0.020, cycles=2):
        state, _ = model.advance_fixed(state, increment, 4)
    samples = np.linspace(0.0, 0.012, 49)
    stresses = []
    for value in samples:
        trial, _ = model.advance_fixed(state, strain_increment(value), 4)
        stresses.append(float(trial.stress[0, 2]))
    differences = np.diff(stresses)
    tolerance = 2.0e-8 * max(1.0, np.max(np.abs(stresses)))
    if np.min(differences) < -tolerance:
        index = int(np.argmin(differences))
        raise AssertionError(
            "nonmonotone local stress map at "
            f"gamma={samples[index]:.6g}: delta_tau={differences[index]:.6g}"
        )


def substep_refinement() -> tuple[float, float]:
    # Stay below the near-liquefaction cone-return regime.  The full
    # stress-controlled qualification below checks large-strain triggering;
    # this local gate isolates constitutive integration error.
    histories = sinusoidal_increments(0.005, points=32, cycles=3)
    states = []
    models = []
    for nsub in (2, 4, 8):
        model = RIVASandMappingBackstressModel()
        models.append(model)
        states.append(initialize(model, void_ratio=0.660, bias_shear=10.0))
    maxima = [0.0, 0.0]
    scale = 1.0
    for increment in histories:
        for index, nsub in enumerate((2, 4, 8)):
            states[index], _ = models[index].advance_fixed(
                states[index], increment, nsub
            )
        scale = max(scale, float(np.linalg.norm(states[2].stress)))
        maxima[0] = max(
            maxima[0], float(np.linalg.norm(states[0].stress - states[1].stress))
        )
        maxima[1] = max(
            maxima[1], float(np.linalg.norm(states[1].stress - states[2].stress))
        )
    errors = maxima[0] / scale, maxima[1] / scale
    if errors[0] > 0.12 or errors[1] > 0.06:
        raise AssertionError(
            f"substep stress errors 2/4={errors[0]:.3%}, 4/8={errors[1]:.3%}"
        )
    return errors


def main() -> None:
    disabled_and_preserved_window_equivalence()
    finite_restart_compatibility_and_outer_surface()
    monotone_local_stress_map()
    errors = substep_refinement()
    print(
        "PASS: preserved-window equivalence, finite/restart/compatibility, "
        "outer-surface correction, monotone stress map, and refinement "
        f"2/4={errors[0]:.3%}, 4/8={errors[1]:.3%}"
    )


if __name__ == "__main__":
    main()
