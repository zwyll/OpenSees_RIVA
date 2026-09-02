"""Intermediate-density biased-flow successor to the mapping checkpoint.

The qualified mapping/backstress mechanism is intentionally retained only in
its original low-density, high-static-bias window.  This research successor
adds a separate reversal-anchored directional compliance law centered on the
PRJ-3484 intermediate-density sand.  The law changes neither the stress
surface nor the compatible volumetric update; it balances plastic flow on the
two biased half-cycles through the shear modulus used inside the local stress
update.

The successor also centers the compatible reversible phase-transformation
wave about a fixed high-bias potential, admits that wave continuously before
the first detected reversal, and ramps it by completed half-cycles.  These
volumetric changes are confined to the same intermediate-density high-bias
window and pressure is always rebuilt from compatible volumetric strain.

The added gates are zero for zero-bias, loose, and dense states.  Setting
``intermediate_bias_flow_enabled=False`` recovers the mapping/backstress
checkpoint exactly.  This is private research code, not production
RIVA-Sand.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from rivasand_port.model import RIVASandState, Tensor, invariants, tensor_norm

from RIVASandMappingBackstress import (
    RIVASandMappingBackstressModel,
    RIVASandMappingBackstressParameters,
    RIVASandMappingBackstressState,
)


@dataclass(frozen=True)
class RIVASandIntermediateBiasFlowParameters(
    RIVASandMappingBackstressParameters
):
    """Mapping checkpoint plus a compact intermediate-density flow law."""

    intermediate_bias_flow_enabled: bool = True
    intermediate_bias_density_onset: float = 0.60
    intermediate_bias_density_peak: float = 0.663
    intermediate_bias_density_cutoff: float = 0.76
    intermediate_bias_onset: float = 0.30
    intermediate_bias_full: float = 0.50
    intermediate_bias_cutoff_onset: float = 0.60
    intermediate_bias_cutoff: float = 0.66

    intermediate_compliance_peak: float = 1.0
    intermediate_compliance_bell_gain: float = 0.50
    intermediate_directional_balance: float = 0.0
    intermediate_balance_bias_exponent: float = 1.0
    intermediate_compliance_minimum: float = 0.20

    intermediate_high_bias_onset: float = 0.62
    intermediate_high_bias_full: float = 0.67
    intermediate_high_bias_cutoff_onset: float = 0.72
    intermediate_high_bias_cutoff: float = 0.78
    intermediate_high_bias_compliance_peak: float = 3.0
    intermediate_high_bias_compliance_bell_gain: float = 1.50
    intermediate_high_bias_directional_balance: float = 0.12
    intermediate_high_bias_amplitude_onset: float = 0.40
    intermediate_high_bias_amplitude_full: float = 0.45
    intermediate_first_reversal_amplitude_ratio: float = 0.40
    intermediate_high_bias_phase_anchor_fraction: float = 0.90
    intermediate_high_bias_phase_relaxation_multiplier: float = 0.25
    intermediate_high_bias_phase_activation_reversals: float = 6.0
    intermediate_high_bias_pre_reversal_phase_scale: float = 1.50
    bias_reversible_volume_ep_ref: float = 1.5e-5

    def __post_init__(self) -> None:
        super().__post_init__()
        if not (
            0.0 <= self.intermediate_bias_density_onset
            < self.intermediate_bias_density_peak
            < self.intermediate_bias_density_cutoff <= 1.0
        ):
            raise ValueError("invalid intermediate-density flow interval")
        if not 0.0 <= self.intermediate_bias_onset < self.intermediate_bias_full:
            raise ValueError("invalid intermediate-bias flow interval")
        if not (
            self.intermediate_bias_full
            < self.intermediate_bias_cutoff_onset
            < self.intermediate_bias_cutoff
        ):
            raise ValueError("invalid intermediate-bias cutoff interval")
        if min(
            self.intermediate_compliance_peak,
            self.intermediate_compliance_bell_gain,
            self.intermediate_high_bias_compliance_peak,
            self.intermediate_high_bias_compliance_bell_gain,
        ) < 0.0:
            raise ValueError("intermediate compliance controls must be nonnegative")
        if not (
            np.isfinite(self.intermediate_directional_balance)
            and np.isfinite(self.intermediate_high_bias_directional_balance)
        ):
            raise ValueError("intermediate directional balance must be finite")
        if self.intermediate_balance_bias_exponent < 0.0:
            raise ValueError("intermediate balance exponent must be nonnegative")
        if not 0.0 < self.intermediate_compliance_minimum <= 1.0:
            raise ValueError("intermediate compliance minimum must lie in (0,1]")
        if not (
            0.0 <= self.intermediate_high_bias_onset
            < self.intermediate_high_bias_full
            < self.intermediate_high_bias_cutoff_onset
            < self.intermediate_high_bias_cutoff
        ):
            raise ValueError("invalid intermediate high-bias interval")
        if not (
            0.0 <= self.intermediate_high_bias_amplitude_onset
            < self.intermediate_high_bias_amplitude_full
        ):
            raise ValueError("invalid intermediate high-bias amplitude interval")
        if not 0.0 < self.intermediate_first_reversal_amplitude_ratio <= 1.0:
            raise ValueError("first-reversal amplitude ratio must lie in (0,1]")
        if not 0.0 <= self.intermediate_high_bias_phase_anchor_fraction <= 2.0:
            raise ValueError("high-bias phase anchor fraction must lie in [0,2]")
        if self.intermediate_high_bias_phase_relaxation_multiplier < 0.0:
            raise ValueError("high-bias phase relaxation multiplier must be nonnegative")
        if self.intermediate_high_bias_phase_activation_reversals <= 0.0:
            raise ValueError("high-bias phase activation reversals must be positive")
        if self.intermediate_high_bias_pre_reversal_phase_scale < 0.0:
            raise ValueError("pre-reversal phase scale must be nonnegative")
        if self.bias_reversible_volume_ep_ref <= 0.0:
            raise ValueError("plastic-activity gate reference must be positive")


@dataclass
class RIVASandIntermediateBiasFlowState(RIVASandMappingBackstressState):
    """Caches immutable density--bias gates at cyclic activation."""

    initial_relative_density_value: float = -1.0
    intermediate_low_gate_value: float = 0.0
    intermediate_high_gate_base: float = 0.0
    ep_half_last: float = 0.0

    def copy(self) -> "RIVASandIntermediateBiasFlowState":
        values: dict[str, object] = {}
        for item in fields(self):
            value = getattr(self, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        return RIVASandIntermediateBiasFlowState(**values)


class RIVASandIntermediateBiasFlowModel(RIVASandMappingBackstressModel):
    """Mapping checkpoint with isolated intermediate biased branch balance."""

    parameters: RIVASandIntermediateBiasFlowParameters
    state: RIVASandIntermediateBiasFlowState | None

    def __init__(
        self, parameters: RIVASandIntermediateBiasFlowParameters | None = None
    ):
        super().__init__(parameters or RIVASandIntermediateBiasFlowParameters())
        # Parameter-only reference density is immutable for the material.
        self._reference_relative_density_value = float(
            super().reference_relative_density()
        )

    def initialize(
        self, stress: Tensor, void_ratio: float
    ) -> RIVASandIntermediateBiasFlowState:
        base = super().initialize(stress, void_ratio)
        values: dict[str, object] = {}
        for item in fields(RIVASandMappingBackstressState):
            value = getattr(base, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        self.state = RIVASandIntermediateBiasFlowState(**values)
        # Initial density depends only on the admitted pressure anchor and the
        # initial state, neither of which evolves during cyclic integration.
        self.state.initial_relative_density_value = float(
            super().initial_relative_density(self.state)
        )
        return self.state.copy()

    def _backbone_forward_euler(
        self,
        old: RIVASandIntermediateBiasFlowState,
        deps: Tensor,
        *,
        force_reversal: bool,
        allow_legacy_reversal: bool,
    ) -> RIVASandIntermediateBiasFlowState:
        state = super()._backbone_forward_euler(
            old,
            deps,
            force_reversal=force_reversal,
            allow_legacy_reversal=allow_legacy_reversal,
        )
        if state.reversals > old.reversals:
            state.ep_half_last = float(old.ep_eq_since_reversal)
        return state

    def inherited_bias_reversible_volume_target(
        self, state: RIVASandState
    ) -> float:
        target = super().inherited_bias_reversible_volume_target(state)
        ep_half_last = float(getattr(state, "ep_half_last", 0.0))
        gate = self._smoothstep(
            ep_half_last / self.parameters.bias_reversible_volume_ep_ref
        )
        return float(target * gate)

    def initial_relative_density(self, state: RIVASandState) -> float:
        if isinstance(state, RIVASandIntermediateBiasFlowState):
            cached = state.initial_relative_density_value
            if cached >= 0.0:
                return float(cached)
        return float(super().initial_relative_density(state))

    def reference_relative_density(self) -> float:
        cached = getattr(self, "_reference_relative_density_value", None)
        if cached is not None:
            return float(cached)
        return float(super().reference_relative_density())

    def _intermediate_density_gate(self, state: RIVASandState) -> float:
        cfg = self.parameters
        density = self.initial_relative_density(state)
        if density <= cfg.intermediate_bias_density_peak:
            return self._smoothstep(
                (density - cfg.intermediate_bias_density_onset)
                / (
                    cfg.intermediate_bias_density_peak
                    - cfg.intermediate_bias_density_onset
                )
            )
        return float(1.0 - self._smoothstep(
            (density - cfg.intermediate_bias_density_peak)
            / (
                cfg.intermediate_bias_density_cutoff
                - cfg.intermediate_bias_density_peak
            )
        ))

    def _raw_intermediate_gate_values(
        self, state: RIVASandState
    ) -> tuple[float, float]:
        cfg = self.parameters
        if not cfg.intermediate_bias_flow_enabled:
            return 0.0, 0.0
        density_gate = self._intermediate_density_gate(state)
        low_bias_gate = self._smoothstep(
            (state.static_bias_index - cfg.intermediate_bias_onset)
            / (cfg.intermediate_bias_full - cfg.intermediate_bias_onset)
        )
        low_bias_gate *= 1.0 - self._smoothstep(
            (state.static_bias_index - cfg.intermediate_bias_cutoff_onset)
            / (
                cfg.intermediate_bias_cutoff
                - cfg.intermediate_bias_cutoff_onset
            )
        )
        high_bias_gate = self._smoothstep(
            (state.static_bias_index - cfg.intermediate_high_bias_onset)
            / (
                cfg.intermediate_high_bias_full
                - cfg.intermediate_high_bias_onset
            )
        )
        high_bias_gate *= 1.0 - self._smoothstep(
            (state.static_bias_index - cfg.intermediate_high_bias_cutoff_onset)
            / (
                cfg.intermediate_high_bias_cutoff
                - cfg.intermediate_high_bias_cutoff_onset
            )
        )
        return float(density_gate * low_bias_gate), float(
            density_gate * high_bias_gate
        )

    def begin_cyclic_phase(
        self, *, reference_stress: Tensor | None = None
    ) -> RIVASandIntermediateBiasFlowState:
        super().begin_cyclic_phase(reference_stress=reference_stress)
        cfg = self.parameters
        if not isinstance(self.state, RIVASandIntermediateBiasFlowState):
            raise TypeError("intermediate biased-flow activation lost its state type")
        low, high = self._raw_intermediate_gate_values(self.state)
        self.state.intermediate_low_gate_value = low
        self.state.intermediate_high_gate_base = high
        # The parent anchors the reversible PT wave at the admitted static
        # potential.  At high bias this makes the wave almost one-sided.  A
        # lower, fixed anchor admits continuous dilation and contraction about
        # the cyclic path while pressure still follows compatible volume.
        anchor_ratio = (
            1.0
            + high
            * (cfg.intermediate_high_bias_phase_anchor_fraction - 1.0)
        )
        self.state.phase_potential_anchor = float(
            self.state.phase_potential_anchor * anchor_ratio
        )
        return self.state.copy()

    def intermediate_bias_flow_gate(self, state: RIVASandState) -> float:
        cfg = self.parameters
        if (
            not cfg.intermediate_bias_flow_enabled
            or not state.cyclic_phase_active
        ):
            return 0.0
        if isinstance(state, RIVASandIntermediateBiasFlowState):
            return float(state.intermediate_low_gate_value)
        return self._raw_intermediate_gate_values(state)[0]

    def intermediate_high_bias_flow_gate(self, state: RIVASandState) -> float:
        cfg = self.parameters
        if (
            not cfg.intermediate_bias_flow_enabled
            or not state.cyclic_phase_active
        ):
            return 0.0
        bias_gate = (
            state.intermediate_high_gate_base
            if isinstance(state, RIVASandIntermediateBiasFlowState)
            else self._raw_intermediate_gate_values(state)[1]
        )
        if bias_gate == 0.0:
            return 0.0
        # The common amplitude memory is a one-sided excursion at the first
        # reversal and a full peak-to-peak range thereafter.  Put both on the
        # same full-cycle scale before this high-bias gate is evaluated.  This
        # prevents the first half-cycle from spuriously activating a mechanism
        # that the stabilized cyclic amplitude would leave inactive.
        stable_amplitude = float(state.cyclic_amplitude)
        if state.amplitude_reversals == 1:
            stable_amplitude *= cfg.intermediate_first_reversal_amplitude_ratio
        amplitude_gate = self._smoothstep(
            (
                stable_amplitude
                - cfg.intermediate_high_bias_amplitude_onset
            )
            / (
                cfg.intermediate_high_bias_amplitude_full
                - cfg.intermediate_high_bias_amplitude_onset
            )
        )
        return float(bias_gate * amplitude_gate)

    def intermediate_branch_multiplier(self, state: RIVASandState) -> float:
        cfg = self.parameters
        low_gate = self.intermediate_bias_flow_gate(state)
        high_gate = self.intermediate_high_bias_flow_gate(state)
        if max(low_gate, high_gate) <= 1.0e-14:
            return 1.0
        transition = self.transformation_progress(state)
        if transition <= 1.0e-14:
            return 1.0

        direction_norm = tensor_norm(state.cyclic_direction)
        projection = float(
            np.sum(state.static_bias_tensor * state.cyclic_direction)
            / max(direction_norm, 1.0e-14)
        )
        signed_bias = float(
            np.sign(projection)
            * (
                abs(projection)
                / max(cfg.branch_compliance_bias_reference, 1.0e-14)
            ) ** cfg.intermediate_balance_bias_exponent
        )
        low_balance = np.exp(
            cfg.intermediate_directional_balance * signed_bias
        )
        high_balance = np.exp(
            cfg.intermediate_high_bias_directional_balance * signed_bias
        )
        zone = self.transformation_zone(state)
        compliance = low_gate * (
            cfg.intermediate_compliance_peak * transition * low_balance
            + cfg.intermediate_compliance_bell_gain
            * zone
        ) + high_gate * (
            cfg.intermediate_high_bias_compliance_peak
            * transition
            * high_balance
            + cfg.intermediate_high_bias_compliance_bell_gain * zone
        )
        return float(np.clip(
            1.0 / (1.0 + compliance),
            cfg.intermediate_compliance_minimum,
            1.0,
        ))

    def phase_reversible_relaxation_multiplier(
        self, state: RIVASandState
    ) -> float:
        gate = self.intermediate_high_bias_flow_gate(state)
        return float(
            1.0
            + gate
            * (
                self.parameters.intermediate_high_bias_phase_relaxation_multiplier
                - 1.0
            )
        )

    def intermediate_high_bias_phase_activation(
        self, state: RIVASandState
    ) -> float:
        if isinstance(state, RIVASandIntermediateBiasFlowState):
            gate = state.intermediate_high_gate_base
        else:
            gate = self._raw_intermediate_gate_values(state)[1]
        if gate == 0.0:
            return 1.0
        ramp = self._smoothstep(
            state.amplitude_reversals
            / self.parameters.intermediate_high_bias_phase_activation_reversals
        )
        return float(1.0 + gate * (ramp - 1.0))

    def _pre_reversal_cyclic_excursion(self, state: RIVASandState) -> float:
        deviator, _, _ = invariants(state.stress)
        static_deviator = state.geostatic_deviator + (
            state.pressure_anchor * state.static_bias_tensor
        )
        return float(
            tensor_norm(deviator - static_deviator)
            / max(state.pressure_anchor, self.parameters.p_min)
        )

    def phase_volume_gate(self, state: RIVASandState) -> float:
        gate = super().phase_volume_gate(state)
        if state.cyclic_phase_active and state.amplitude_reversals < 1:
            high_gate = (
                state.intermediate_high_gate_base
                if isinstance(state, RIVASandIntermediateBiasFlowState)
                else self._raw_intermediate_gate_values(state)[1]
            )
            gate = max(
                gate,
                high_gate * self.phase_volume_density_weight(state),
            )
        return float(gate)

    def phase_activity(
        self,
        state: RIVASandState,
        *,
        bias_exponent: float | None = None,
        pressure_exponent: float | None = None,
    ) -> float:
        activity = super().phase_activity(
            state,
            bias_exponent=bias_exponent,
            pressure_exponent=pressure_exponent,
        )
        if activity <= 0.0 and state.cyclic_phase_active:
            high_gate = (
                state.intermediate_high_gate_base
                if isinstance(state, RIVASandIntermediateBiasFlowState)
                else self._raw_intermediate_gate_values(state)[1]
            )
            if high_gate > 1.0e-14 and state.amplitude_reversals < 1:
                cfg = self.parameters
                _, pressure, _ = invariants(state.stress)
                excursion = self._pre_reversal_cyclic_excursion(state)
                ratio = excursion / cfg.cyclic_amplitude_reference
                amplitude_activity = self._smoothstep(
                    (ratio - cfg.phase_amplitude_onset_ratio)
                    / (
                        cfg.phase_amplitude_full_ratio
                        - cfg.phase_amplitude_onset_ratio
                    )
                )
                exponent = (
                    cfg.phase_pressure_exponent
                    if pressure_exponent is None
                    else pressure_exponent
                )
                pressure_ratio = np.clip(
                    pressure / max(state.pressure_anchor, cfg.p_min),
                    0.0,
                    1.0,
                ) ** exponent
                bias = max(state.static_bias_index, 0.0)
                reference = cfg.branch_compliance_bias_reference
                onset = 0.25 * reference
                bias_gate = self._smoothstep(
                    bias / max(onset, 1.0e-14)
                )
                power = (
                    cfg.phase_bias_exponent
                    if bias_exponent is None
                    else bias_exponent
                )
                bias_activity = min(
                    bias_gate * (reference / max(bias, 1.0e-14)) ** power,
                    2.0,
                )
                activity = float(
                    high_gate
                    * cfg.intermediate_high_bias_pre_reversal_phase_scale
                    * self.phase_volume_density_weight(state)
                    * bias_activity
                    * amplitude_activity
                    * pressure_ratio
                )
                return activity
        return float(
            activity * self.intermediate_high_bias_phase_activation(state)
        )

    def cyclic_flow_factors(
        self, pressure: float, state: RIVASandState
    ) -> tuple[float, float]:
        """Evaluate the inherited three-level cyclic-flow blend once.

        The base, phase-transformation, and loose-flow implementations each
        recompute the same pressure ratio and power before blending their
        shear factors.  This specialization preserves their operation order
        while sharing that common scalar activity.
        """
        cfg = self.parameters
        if (
            not cfg.cyclic_flow_correction_enabled
            or state.amplitude_reversals < cfg.cyclic_flow_minimum_reversals
        ):
            return 1.0, 1.0
        pressure_ratio = np.clip(
            pressure / max(state.pressure_anchor, cfg.p_min), 0.0, 1.0
        )
        activity = pressure_ratio**cfg.cyclic_flow_pressure_exponent
        shear = 1.0 - cfg.cyclic_shear_modulus_reduction * activity
        hardening = 1.0 + cfg.cyclic_hardening_boost * activity
        if cfg.phase_transformation_enabled:
            phase_shear = (
                1.0
                - cfg.phase_cyclic_shear_modulus_reduction * activity
            )
            phase_gate = self.phase_dense_bias_gate(state)
            shear = shear + phase_gate * (phase_shear - shear)
        if cfg.loose_shear_flow_enabled:
            loose_gate = self.loose_phase_gate(state)
            if loose_gate > 1.0e-14:
                loose_shear = (
                    1.0
                    - cfg.loose_shear_cyclic_modulus_reduction * activity
                )
                shear = shear + loose_gate * (loose_shear - shear)
        return float(shear), float(hardening)

    def moduli_for_state(
        self, pressure: float, state: RIVASandState
    ) -> tuple[float, float]:
        shear, bulk = super().moduli_for_state(pressure, state)
        shear *= self.intermediate_branch_multiplier(state)
        return float(shear), float(bulk)

    def dss_history_values(
        self, state: RIVASandState | None = None
    ) -> dict[str, float]:
        current = state or self.state
        if current is None:
            raise RuntimeError("initialize the model first")
        values = super().dss_history_values(current)
        values.update(
            intermediate_bias_flow_gate=self.intermediate_bias_flow_gate(current),
            intermediate_high_bias_flow_gate=(
                self.intermediate_high_bias_flow_gate(current)
            ),
            intermediate_branch_multiplier=self.intermediate_branch_multiplier(
                current
            ),
            intermediate_stable_amplitude=(
                float(current.cyclic_amplitude)
                * (
                    self.parameters.intermediate_first_reversal_amplitude_ratio
                    if current.amplitude_reversals == 1
                    else 1.0
                )
            ),
            intermediate_phase_potential_anchor=float(
                getattr(current, "phase_potential_anchor", 0.0)
            ),
            intermediate_high_bias_phase_activation=(
                self.intermediate_high_bias_phase_activation(current)
            ),
            intermediate_pre_reversal_cyclic_excursion=(
                self._pre_reversal_cyclic_excursion(current)
            ),
        )
        return values


__all__ = [
    "RIVASandIntermediateBiasFlowModel",
    "RIVASandIntermediateBiasFlowParameters",
    "RIVASandIntermediateBiasFlowState",
]
