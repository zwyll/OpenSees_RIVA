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
    @staticmethod
    def drive_cycles(
        model: RIVASandIntermediateBiasFlowModel,
        amplitude: float,
        cycles: int = 2,
    ) -> None:
        previous = 0.0
        for step in range(1, 32 * cycles + 1):
            gamma = amplitude * np.sin(2.0 * np.pi * step / 32.0)
            deps = np.zeros((3, 3))
            deps[0, 2] = deps[2, 0] = 0.5 * (gamma - previous)
            assert model.state is not None
            model.state, _ = model.advance_fixed(model.state, deps, 4)
            previous = gamma

    def test_initial_relative_density_is_cached_and_copied(self) -> None:
        model = RIVASandIntermediateBiasFlowModel()
        state = model.initialize(initial_stress(10.0), 0.601)
        cached = state.initial_relative_density_value
        self.assertGreaterEqual(cached, 0.0)
        self.assertEqual(model.initial_relative_density(state), cached)
        self.assertEqual(state.copy().initial_relative_density_value, cached)

    def test_disabled_intermediate_branch_matches_preserved_mapping_path(self) -> None:
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

    def test_activation_does_not_fabricate_reversal_history(self) -> None:
        model = RIVASandIntermediateBiasFlowModel()
        model.initialize(initial_stress(10.0), 0.601)
        model.begin_cyclic_phase(reference_stress=initial_stress(0.0))
        assert model.state is not None
        self.assertEqual(model.state.amplitude_reversals, 0)
        self.assertEqual(model.state.ep_half_last, 0.0)

    def test_microcycle_suppresses_inherited_bias_wave(self) -> None:
        model = RIVASandIntermediateBiasFlowModel(replace(
            RIVASandIntermediateBiasFlowParameters(),
            phase_transformation_enabled=False,
            mapping_backstress_enabled=False,
            intermediate_bias_flow_enabled=False,
        ))
        model.initialize(initial_stress(10.0), 0.601)
        model.begin_cyclic_phase(reference_stress=initial_stress(0.0))
        self.drive_cycles(model, 1.0e-5)
        assert model.state is not None
        self.assertLess(
            model.state.ep_half_last,
            model.parameters.bias_reversible_volume_ep_ref,
        )
        self.assertLess(abs(model.state.bias_reversible_volume), 1.0e-10)

    def test_both_backbones_capture_completed_plastic_activity(self) -> None:
        reference = initial_stress(0.0)
        for void_ratio, shear, mapping_active in (
            (0.601, 10.0, False),
            (0.662, 15.0, True),
        ):
            with self.subTest(mapping_active=mapping_active):
                model = RIVASandIntermediateBiasFlowModel()
                model.initialize(initial_stress(shear), void_ratio)
                model.begin_cyclic_phase(reference_stress=reference)
                assert model.state is not None
                self.assertEqual(
                    model.mapping_gate(model.state) > 1.0e-14,
                    mapping_active,
                )
                self.drive_cycles(model, 0.01)
                self.assertGreater(
                    model.state.ep_half_last,
                    model.parameters.bias_reversible_volume_ep_ref,
                )

    def test_plastic_activity_refines_from_four_to_eight_substeps(self) -> None:
        reference = initial_stress(0.0)
        for void_ratio, shear in ((0.601, 10.0), (0.662, 15.0)):
            with self.subTest(void_ratio=void_ratio, shear=shear):
                histories: list[np.ndarray] = []
                states = []
                for nsub in (4, 8):
                    model = RIVASandIntermediateBiasFlowModel()
                    model.initialize(initial_stress(shear), void_ratio)
                    model.begin_cyclic_phase(reference_stress=reference)
                    previous = 0.0
                    history = []
                    for step in range(1, 129):
                        gamma = 0.01 * np.sin(2.0 * np.pi * step / 32.0)
                        deps = np.zeros((3, 3))
                        deps[0, 2] = deps[2, 0] = 0.5 * (gamma - previous)
                        assert model.state is not None
                        model.state, _ = model.advance_fixed(
                            model.state, deps, nsub
                        )
                        history.append(model.state.stress.copy())
                        previous = gamma
                    histories.append(np.asarray(history))
                    states.append(model.state.copy())
                stress_error = float(np.max(np.linalg.norm(
                    histories[0] - histories[1], axis=(1, 2)
                )) / np.max(np.linalg.norm(histories[1], axis=(1, 2))))
                ep_error = (
                    abs(states[0].ep_half_last - states[1].ep_half_last)
                    / states[1].ep_half_last
                )
                self.assertLess(stress_error, 0.03)
                self.assertLess(ep_error, 0.02)
                self.assertEqual(
                    states[0].amplitude_reversals,
                    states[1].amplitude_reversals,
                )


if __name__ == "__main__":
    unittest.main()
