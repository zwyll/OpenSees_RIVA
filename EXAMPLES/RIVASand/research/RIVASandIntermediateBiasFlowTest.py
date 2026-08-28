"""Focused regression tests for the intermediate biased-flow successor."""

from __future__ import annotations

from dataclasses import fields, replace
from pathlib import Path
import sys
import unittest

import numpy as np


HERE = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(HERE))

from RIVASandIntermediateBiasFlow import (  # noqa: E402
    RIVASandIntermediateBiasFlowModel,
    RIVASandIntermediateBiasFlowParameters,
)
from RIVASandMappingBackstress import (  # noqa: E402
    RIVASandMappingBackstressModel,
    RIVASandMappingBackstressParameters,
    RIVASandMappingBackstressState,
)


def initial_stress(shear: float) -> np.ndarray:
    return np.asarray([
        [-19.4, 0.0, shear],
        [0.0, -19.4, 0.0],
        [shear, 0.0, -40.0],
    ])


class IntermediateBiasFlowTest(unittest.TestCase):
    def test_initial_relative_density_is_cached_and_copied(self) -> None:
        model = RIVASandIntermediateBiasFlowModel()
        state = model.initialize(initial_stress(10.0), 0.601)
        cached = state.initial_relative_density_value
        self.assertGreaterEqual(cached, 0.0)
        self.assertEqual(model.initial_relative_density(state), cached)
        self.assertEqual(state.copy().initial_relative_density_value, cached)

    def test_disabled_successor_is_exact_mapping_checkpoint(self) -> None:
        parent = RIVASandMappingBackstressModel(
            RIVASandMappingBackstressParameters()
        )
        child = RIVASandIntermediateBiasFlowModel(replace(
            RIVASandIntermediateBiasFlowParameters(),
            intermediate_bias_flow_enabled=False,
        ))
        reference = initial_stress(0.0)
        for model in (parent, child):
            model.initialize(initial_stress(10.0), 0.601)
            model.begin_cyclic_phase(reference_stress=reference)
        increments = (0.0004, -0.0008, 0.0008, -0.0008, 0.0004)
        for engineering_shear in increments:
            deps = np.zeros((3, 3))
            deps[0, 2] = deps[2, 0] = engineering_shear / 2.0
            parent.state, _ = parent.advance_fixed(parent.state, deps, 4)
            child.state, _ = child.advance_fixed(child.state, deps, 4)
        assert parent.state is not None and child.state is not None
        for item in fields(RIVASandMappingBackstressState):
            expected = getattr(parent.state, item.name)
            actual = getattr(child.state, item.name)
            if isinstance(expected, np.ndarray):
                self.assertTrue(np.array_equal(actual, expected), item.name)
            else:
                self.assertEqual(actual, expected, item.name)

    def test_bias_windows_are_nonoverlapping(self) -> None:
        model = RIVASandIntermediateBiasFlowModel()
        reference = initial_stress(0.0)

        model.initialize(initial_stress(10.35), 0.601)
        model.begin_cyclic_phase(reference_stress=reference)
        assert model.state is not None
        model.state.cyclic_amplitude = 0.65
        model.state.amplitude_reversals = 2
        self.assertGreater(model.intermediate_bias_flow_gate(model.state), 0.95)
        self.assertLess(model.intermediate_high_bias_flow_gate(model.state), 1.0e-12)

        model.initialize(initial_stress(12.60), 0.601)
        model.begin_cyclic_phase(reference_stress=reference)
        assert model.state is not None
        model.state.cyclic_amplitude = 0.65
        model.state.amplitude_reversals = 2
        self.assertLess(model.intermediate_bias_flow_gate(model.state), 1.0e-12)
        self.assertGreater(model.intermediate_high_bias_flow_gate(model.state), 0.95)

    def test_first_reversal_uses_full_cycle_amplitude_scale(self) -> None:
        model = RIVASandIntermediateBiasFlowModel()
        reference = initial_stress(0.0)
        model.initialize(initial_stress(12.6), 0.601)
        model.begin_cyclic_phase(reference_stress=reference)
        assert model.state is not None
        model.state.cyclic_amplitude = 0.90
        model.state.amplitude_reversals = 1
        self.assertEqual(model.intermediate_high_bias_flow_gate(model.state), 0.0)
        model.state.cyclic_amplitude = 0.65
        model.state.amplitude_reversals = 2
        self.assertGreater(model.intermediate_high_bias_flow_gate(model.state), 0.95)

    def test_high_bias_phase_volume_is_admitted_before_first_reversal(self) -> None:
        model = RIVASandIntermediateBiasFlowModel()
        reference = initial_stress(0.0)
        model.initialize(initial_stress(12.6), 0.601)
        model.begin_cyclic_phase(reference_stress=reference)
        assert model.state is not None
        self.assertEqual(model.state.amplitude_reversals, 0)
        self.assertGreater(model.phase_volume_gate(model.state), 0.25)

    def test_zero_bias_and_dense_states_are_unchanged(self) -> None:
        model = RIVASandIntermediateBiasFlowModel()
        reference = initial_stress(0.0)
        for shear, void_ratio in ((0.0, 0.601), (10.0, 0.545)):
            model.initialize(initial_stress(shear), void_ratio)
            model.begin_cyclic_phase(reference_stress=reference)
            assert model.state is not None
            model.state.cyclic_amplitude = 0.65
            model.state.amplitude_reversals = 2
            self.assertEqual(model.intermediate_bias_flow_gate(model.state), 0.0)
            self.assertEqual(
                model.intermediate_high_bias_flow_gate(model.state), 0.0
            )
            self.assertEqual(model.intermediate_branch_multiplier(model.state), 1.0)


if __name__ == "__main__":
    unittest.main()
