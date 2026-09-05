"""Loose-sand stabilization and zero-bias PT research successor.

This private prototype starts from the restrained accumulation-control
checkpoint and adds two compact, non-overlapping density/bias windows:

* low-density biased stabilization plus a signed phase-volume law for the
  PRJ-4666 loose test; and
* a fixed-direction, symmetric phase-volume law for dense zero-bias loading.

Both additions vanish exactly outside their declared density/bias windows.
This is research code, not production RIVA-Sand.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields

import numpy as np

from rivasand_port.model import RIVASandState, SQRT23, Tensor, invariants, tensor_norm

from RIVASandAccumulationControl import (
    RIVASandAccumulationControlModel,
    RIVASandAccumulationControlParameters,
    RIVASandAccumulationControlState,
)


@dataclass(frozen=True)
class RIVASandLooseUnbiasedCorrectionParameters(
    RIVASandAccumulationControlParameters
):
    """Restrained accumulation plus two isolated DSS correction windows."""

    loose_stabilization_enabled: bool = True
    loose_stabilization_density_full: float = 0.50
    loose_stabilization_density_cutoff: float = 0.56
    loose_stabilization_bias_onset: float = 0.13
    loose_stabilization_hardening_multiplier: float = 1.40
    loose_stabilization_contraction_scale: float = 0.50
    loose_phase_enabled: bool = True
    loose_phase_activity_scale: float = 0.20
    loose_phase_wave_activity_scale: float = 0.30
    loose_phase_mean_activity_scale: float = 1.10
    loose_phase_replacement_fraction: float = 1.0

    unbiased_phase_enabled: bool = True
    unbiased_phase_density_onset: float = 0.82
    unbiased_phase_density_full: float = 0.89
    unbiased_phase_bias_cutoff: float = 0.08
    unbiased_phase_activity_scale: float = 0.70
    unbiased_phase_potential_center: float = 0.80

    def __post_init__(self) -> None:
        super().__post_init__()
        if not (
            0.0 <= self.loose_stabilization_density_full
            < self.loose_stabilization_density_cutoff <= 1.0
        ):
            raise ValueError("invalid loose-stabilization density interval")
        if self.loose_stabilization_bias_onset <= 0.0:
            raise ValueError("loose-stabilization bias onset must be positive")
        if self.loose_stabilization_hardening_multiplier <= 0.0:
            raise ValueError("loose-stabilization hardening multiplier must be positive")
        if not 0.0 <= self.loose_stabilization_contraction_scale <= 1.0:
            raise ValueError("loose-stabilization contraction scale must lie in [0,1]")
        if min(
            self.loose_phase_activity_scale,
            self.loose_phase_wave_activity_scale,
            self.loose_phase_mean_activity_scale,
        ) < 0.0:
            raise ValueError("loose phase activity scales must be nonnegative")
        if not 0.0 <= self.loose_phase_replacement_fraction <= 1.0:
            raise ValueError("loose phase replacement fraction must lie in [0,1]")
        if not (
            0.0 <= self.unbiased_phase_density_onset
            < self.unbiased_phase_density_full <= 1.0
        ):
            raise ValueError("invalid unbiased-phase density interval")
        if self.unbiased_phase_bias_cutoff <= 0.0:
            raise ValueError("unbiased-phase bias cutoff must be positive")
        if self.unbiased_phase_activity_scale < 0.0:
            raise ValueError("unbiased-phase activity scale must be nonnegative")
        if not 0.0 <= self.unbiased_phase_potential_center <= 1.0:
            raise ValueError("unbiased-phase potential center must lie in [0,1]")


@dataclass
class RIVASandLooseUnbiasedCorrectionState(RIVASandAccumulationControlState):
    """Adds one fixed zero-bias phase direction."""

    unbiased_phase_direction: Tensor = field(
        default_factory=lambda: np.zeros((3, 3), dtype=float)
    )

    def copy(self) -> "RIVASandLooseUnbiasedCorrectionState":
        values = {}
        for item in fields(RIVASandLooseUnbiasedCorrectionState):
            value = getattr(self, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        return RIVASandLooseUnbiasedCorrectionState(**values)


class RIVASandLooseUnbiasedCorrectionModel(RIVASandAccumulationControlModel):
    """Latest research checkpoint with isolated loose and zero-bias repairs."""

    parameters: RIVASandLooseUnbiasedCorrectionParameters
    state: RIVASandLooseUnbiasedCorrectionState | None

    def __init__(
        self, parameters: RIVASandLooseUnbiasedCorrectionParameters | None = None
    ):
        super().__init__(parameters or RIVASandLooseUnbiasedCorrectionParameters())

    def initialize(
        self, stress: Tensor, void_ratio: float
    ) -> RIVASandLooseUnbiasedCorrectionState:
        base = super().initialize(stress, void_ratio)
        values = {
            item.name: getattr(base, item.name)
            for item in fields(RIVASandAccumulationControlState)
        }
        for name, value in tuple(values.items()):
            if isinstance(value, np.ndarray):
                values[name] = value.copy()
        self.state = RIVASandLooseUnbiasedCorrectionState(**values)
        return self.state.copy()

    def begin_cyclic_phase(
        self, *, reference_stress: Tensor | None = None
    ) -> RIVASandLooseUnbiasedCorrectionState:
        super().begin_cyclic_phase(reference_stress=reference_stress)
        if not isinstance(self.state, RIVASandLooseUnbiasedCorrectionState):
            raise TypeError("loose/unbiased update lost its state type")
        self.state.unbiased_phase_direction.fill(0.0)
        return self.state.copy()

    def loose_stabilization_gate(self, state: RIVASandState) -> float:
        cfg = self.parameters
        if not cfg.loose_stabilization_enabled:
            return 0.0
        density_position = (
            self.initial_relative_density(state)
            - cfg.loose_stabilization_density_full
        ) / (
            cfg.loose_stabilization_density_cutoff
            - cfg.loose_stabilization_density_full
        )
        density_gate = 1.0 - self._smoothstep(density_position)
        bias_gate = self._smoothstep(
            state.static_bias_index / cfg.loose_stabilization_bias_onset
        )
        return float(density_gate * bias_gate)

    def unbiased_phase_gate(self, state: RIVASandState) -> float:
        cfg = self.parameters
        if not cfg.unbiased_phase_enabled:
            return 0.0
        density_position = (
            self.initial_relative_density(state) - cfg.unbiased_phase_density_onset
        ) / (
            cfg.unbiased_phase_density_full - cfg.unbiased_phase_density_onset
        )
        density_gate = self._smoothstep(density_position)
        bias_gate = 1.0 - self._smoothstep(
            state.static_bias_index / cfg.unbiased_phase_bias_cutoff
        )
        return float(density_gate * bias_gate)

    def loose_phase_gate(self, state: RIVASandState) -> float:
        if not self.parameters.loose_phase_enabled:
            return 0.0
        return self.loose_stabilization_gate(state)

    def hardening_prefactor_for_state(
        self, pressure: float, state: RIVASandState
    ) -> float:
        hardening = super().hardening_prefactor_for_state(pressure, state)
        gate = self.loose_stabilization_gate(state)
        multiplier = 1.0 + gate * (
            self.parameters.loose_stabilization_hardening_multiplier - 1.0
        )
        return float(hardening * multiplier)

    def irreversible_contraction_factor(self, state: RIVASandState) -> float:
        """Restrain only low-density biased volumetric accumulation.

        This leaves the shear mapping, plastic modulus, and loop width intact.
        The correction instead scales the irreversible dilatancy component that
        drives compatible pressure loss.
        """
        factor = super().irreversible_contraction_factor(state)
        gate = self.loose_stabilization_gate(state)
        scale = 1.0 - gate * (
            1.0 - self.parameters.loose_stabilization_contraction_scale
        )
        return float(factor * scale)

    def phase_volume_gate(self, state: RIVASandState) -> float:
        loose = self.loose_phase_gate(state)
        loose_activation = max(
            self.parameters.loose_phase_activity_scale,
            self.parameters.loose_phase_wave_activity_scale,
            self.parameters.loose_phase_mean_activity_scale,
            self.parameters.loose_phase_replacement_fraction,
        )
        return float(max(
            super().phase_volume_gate(state),
            self.unbiased_phase_gate(state),
            loose * loose_activation,
        ))

    def phase_volume_replacement_gate(self, state: RIVASandState) -> float:
        return float(max(
            super().phase_volume_replacement_gate(state),
            self.unbiased_phase_gate(state),
            self.loose_phase_gate(state)
            * self.parameters.loose_phase_replacement_fraction,
        ))

    def phase_activity(
        self,
        state: RIVASandState,
        *,
        bias_exponent: float | None = None,
        pressure_exponent: float | None = None,
    ) -> float:
        biased = super().phase_activity(
            state,
            bias_exponent=bias_exponent,
            pressure_exponent=pressure_exponent,
        )
        unbiased_gate = self.unbiased_phase_gate(state)
        loose_gate = self.loose_phase_gate(state)
        if (
            max(unbiased_gate, loose_gate) <= 1.0e-14
            or not state.cyclic_phase_active
        ):
            return biased
        _, pressure, _ = invariants(state.stress)
        exponent = (
            self.parameters.phase_pressure_exponent
            if pressure_exponent is None
            else pressure_exponent
        )
        pressure_ratio = np.clip(
            pressure / max(state.pressure_anchor, self.parameters.p_min),
            0.0,
            1.0,
        ) ** exponent
        unbiased = (
            self.parameters.unbiased_phase_activity_scale
            * unbiased_gate
            * self.phase_amplitude_activity(state)
            * pressure_ratio
        )
        loose_scale = self.parameters.loose_phase_activity_scale
        if bias_exponent == self.parameters.phase_wave_bias_exponent:
            loose_scale = self.parameters.loose_phase_wave_activity_scale
        elif bias_exponent == self.parameters.phase_mean_bias_exponent:
            loose_scale = self.parameters.loose_phase_mean_activity_scale
        loose = (
            loose_scale
            * loose_gate
            * self.phase_amplitude_activity(state)
            * pressure_ratio
        )
        return float(biased + unbiased + loose)

    def signed_phase_potential(self, state: RIVASandState) -> float:
        if self.unbiased_phase_gate(state) <= 1.0e-14:
            return super().signed_phase_potential(state)
        direction = getattr(state, "unbiased_phase_direction", None)
        if direction is None or tensor_norm(direction) <= 1.0e-14:
            return 0.0
        direction = direction / tensor_norm(direction)
        signed_eta = float(np.sum(state.alpha * direction))
        _, pressure, _ = invariants(state.stress)
        _, md, _ = self.surfaces(pressure, state.void_ratio)
        eta_pt = self.parameters.phase_ratio * SQRT23 * md
        # Dense zero-bias DSS is symmetric: both positive and negative shear
        # peaks dilate relative to the low-stress portion of a cycle.  A smooth
        # absolute coordinate therefore supplies one continuous two-sided
        # target rather than reversing the volumetric target every half-cycle.
        smooth_magnitude = np.sqrt(signed_eta * signed_eta + 1.0e-12) - 1.0e-6
        magnitude_potential = float(np.tanh(
            smooth_magnitude
            / max(self.parameters.phase_width * eta_pt, 1.0e-12)
        ))
        return magnitude_potential - self.parameters.unbiased_phase_potential_center

    def phase_coordinates(
        self, state: RIVASandState
    ) -> tuple[float, float, float]:
        if self.unbiased_phase_gate(state) <= 1.0e-14:
            return super().phase_coordinates(state)
        direction = getattr(state, "unbiased_phase_direction", None)
        if direction is None or tensor_norm(direction) <= 1.0e-14:
            return 0.0, 0.0, -1.0
        direction = direction / tensor_norm(direction)
        loading_eta = abs(float(np.sum(state.alpha * direction)))
        _, pressure, _ = invariants(state.stress)
        _, md, _ = self.surfaces(pressure, state.void_ratio)
        eta_pt = self.parameters.phase_ratio * SQRT23 * md
        signed_phase = float(np.tanh(
            (loading_eta - eta_pt)
            / max(self.parameters.phase_width * eta_pt, 1.0e-12)
        ))
        return loading_eta, float(eta_pt), signed_phase

    def advance_fixed(
        self,
        initial: RIVASandLooseUnbiasedCorrectionState,
        deps: Tensor,
        nsub: int = 1,
    ):
        state, info = super().advance_fixed(initial, deps, nsub)
        if not isinstance(state, RIVASandLooseUnbiasedCorrectionState):
            raise TypeError("loose/unbiased update lost its state type")
        if (
            self.unbiased_phase_gate(state) > 1.0e-14
            and tensor_norm(state.unbiased_phase_direction) <= 1.0e-14
            and tensor_norm(state.cyclic_direction) > 1.0e-14
        ):
            state.unbiased_phase_direction = (
                state.cyclic_direction / tensor_norm(state.cyclic_direction)
            )
        return state, info

    def dss_history_values(
        self, state: RIVASandState | None = None
    ) -> dict[str, float]:
        current = state or self.state
        if current is None:
            raise RuntimeError("initialize the model first")
        values = super().dss_history_values(current)
        values.update(
            loose_stabilization_gate=self.loose_stabilization_gate(current),
            loose_phase_gate=self.loose_phase_gate(current),
            unbiased_phase_gate=self.unbiased_phase_gate(current),
            unbiased_phase_direction_norm=tensor_norm(
                getattr(current, "unbiased_phase_direction", np.zeros((3, 3)))
            ),
        )
        return values


__all__ = [
    "RIVASandLooseUnbiasedCorrectionModel",
    "RIVASandLooseUnbiasedCorrectionParameters",
    "RIVASandLooseUnbiasedCorrectionState",
]
