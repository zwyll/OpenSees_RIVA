"""Plastic-work shakedown successor to the restrained PT research checkpoint.

This private prototype adds one intermediate-density hardening memory to the
restrained phase-transformation model.  The memory is the accumulated plastic
multiplier after cyclic activation, not elapsed time or an external cycle
counter.  Its amplitude-dependent onset preserves the calibrated early loop
at moderate CSR while activating earlier in the high-CSR biased test.

Setting ``phase_accumulation_control_enabled=False`` recovers the restrained
phase-transformation checkpoint exactly.  This is not production RIVA-Sand.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from rivasand_port.model import RIVASandState, Tensor

from RIVASandPhaseTransformation import (
    RIVASandPhaseTransformationModel,
    RIVASandPhaseTransformationParameters,
    RIVASandPhaseTransformationState,
)


@dataclass(frozen=True)
class RIVASandAccumulationControlParameters(
    RIVASandPhaseTransformationParameters
):
    """Restrained PT calibration plus smooth plastic-work shakedown."""

    phase_accumulation_control_enabled: bool = True
    phase_accumulation_density_onset: float = 0.53
    phase_accumulation_density_peak: float = 0.663
    phase_accumulation_density_cutoff: float = 0.82
    phase_accumulation_reference_amplitude: float = 0.43
    phase_accumulation_memory_onset: float = 0.030
    phase_accumulation_memory_width: float = 0.006
    phase_accumulation_amplitude_exponent: float = 1.0
    phase_accumulation_gain_full_amplitude: float = 0.65
    phase_accumulation_high_amplitude_gain_ratio: float = 0.6666666666666666
    phase_accumulation_hardening_gain: float = 3.0

    def __post_init__(self) -> None:
        super().__post_init__()
        if not (
            0.0 <= self.phase_accumulation_density_onset
            < self.phase_accumulation_density_peak
            < self.phase_accumulation_density_cutoff <= 1.0
        ):
            raise ValueError("invalid accumulation-control density interval")
        if self.phase_accumulation_reference_amplitude <= 0.0:
            raise ValueError("accumulation reference amplitude must be positive")
        if self.phase_accumulation_memory_onset <= 0.0:
            raise ValueError("accumulation memory onset must be positive")
        if self.phase_accumulation_memory_width <= 0.0:
            raise ValueError("accumulation memory width must be positive")
        if self.phase_accumulation_amplitude_exponent < 0.0:
            raise ValueError("accumulation amplitude exponent must be nonnegative")
        if (
            self.phase_accumulation_gain_full_amplitude
            <= self.phase_accumulation_reference_amplitude
        ):
            raise ValueError("accumulation gain full amplitude must exceed reference")
        if not 0.0 < self.phase_accumulation_high_amplitude_gain_ratio <= 1.0:
            raise ValueError("high-amplitude accumulation gain ratio must lie in (0,1]")
        if self.phase_accumulation_hardening_gain < 0.0:
            raise ValueError("accumulation hardening gain must be nonnegative")


@dataclass
class RIVASandAccumulationControlState(RIVASandPhaseTransformationState):
    """Adds cyclic-work memory and its host-step committed hardening state."""

    phase_accumulation_lambda_anchor: float = 0.0
    phase_accumulation_hardening_state: float = 0.0

    def copy(self) -> "RIVASandAccumulationControlState":
        values = {}
        for item in fields(RIVASandAccumulationControlState):
            value = getattr(self, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        return RIVASandAccumulationControlState(**values)


class RIVASandAccumulationControlModel(RIVASandPhaseTransformationModel):
    """Restrained PT model with intermediate-density accumulation control."""

    parameters: RIVASandAccumulationControlParameters
    state: RIVASandAccumulationControlState | None

    def __init__(
        self, parameters: RIVASandAccumulationControlParameters | None = None
    ):
        super().__init__(parameters or RIVASandAccumulationControlParameters())

    def initialize(
        self, stress: Tensor, void_ratio: float
    ) -> RIVASandAccumulationControlState:
        base = super().initialize(stress, void_ratio)
        values = {
            item.name: getattr(base, item.name)
            for item in fields(RIVASandPhaseTransformationState)
        }
        for name, value in tuple(values.items()):
            if isinstance(value, np.ndarray):
                values[name] = value.copy()
        self.state = RIVASandAccumulationControlState(**values)
        return self.state.copy()

    def begin_cyclic_phase(
        self, *, reference_stress: Tensor | None = None
    ) -> RIVASandAccumulationControlState:
        super().begin_cyclic_phase(reference_stress=reference_stress)
        if not isinstance(self.state, RIVASandAccumulationControlState):
            raise TypeError("accumulation-control update lost its state type")
        self.state.phase_accumulation_lambda_anchor = float(self.state.lambda_total)
        self.state.phase_accumulation_hardening_state = 0.0
        return self.state.copy()

    def phase_accumulation_density_gate(self, state: RIVASandState) -> float:
        """Compact density window centered on the PRJ-3484 biased sand."""
        cfg = self.parameters
        density = self.initial_relative_density(state)
        if density <= cfg.phase_accumulation_density_peak:
            position = (
                density - cfg.phase_accumulation_density_onset
            ) / (
                cfg.phase_accumulation_density_peak
                - cfg.phase_accumulation_density_onset
            )
            return self._smoothstep(position)
        position = (
            density - cfg.phase_accumulation_density_peak
        ) / (
            cfg.phase_accumulation_density_cutoff
            - cfg.phase_accumulation_density_peak
        )
        return float(1.0 - self._smoothstep(position))

    def phase_accumulation_bias_gate(self, state: RIVASandState) -> float:
        """Make the shakedown memory vanish continuously at zero static bias."""
        bias = self.projected_bias(state)
        onset = 0.50 * self.parameters.branch_compliance_bias_reference
        return self._smoothstep(bias / max(onset, 1.0e-14))

    def phase_accumulation_memory(self, state: RIVASandState) -> float:
        anchor = float(getattr(state, "phase_accumulation_lambda_anchor", 0.0))
        return float(max(state.lambda_total - anchor, 0.0))

    def phase_accumulation_onset(self, state: RIVASandState) -> float:
        """Lower the plastic-work onset smoothly as cyclic amplitude rises."""
        cfg = self.parameters
        amplitude = max(
            state.cyclic_amplitude,
            cfg.phase_accumulation_reference_amplitude,
        )
        ratio = cfg.phase_accumulation_reference_amplitude / amplitude
        return float(
            cfg.phase_accumulation_memory_onset
            * ratio**cfg.phase_accumulation_amplitude_exponent
        )

    def phase_accumulation_target_activity(self, state: RIVASandState) -> float:
        """End-of-host-step activity targeted by the plastic-work memory."""
        cfg = self.parameters
        if (
            not cfg.phase_transformation_enabled
            or not cfg.phase_accumulation_control_enabled
            or not state.cyclic_phase_active
        ):
            return 0.0
        onset = self.phase_accumulation_onset(state)
        transition = self._smoothstep(
            (self.phase_accumulation_memory(state) - onset)
            / cfg.phase_accumulation_memory_width
        )
        return float(
            self.phase_accumulation_density_gate(state)
            * self.phase_accumulation_bias_gate(state)
            * transition
        )

    def phase_accumulation_amplitude_gain(self, state: RIVASandState) -> float:
        """Reduce late-cycle hardening smoothly at the larger CSR amplitude."""
        cfg = self.parameters
        amplitude_position = (
            state.cyclic_amplitude - cfg.phase_accumulation_reference_amplitude
        ) / (
            cfg.phase_accumulation_gain_full_amplitude
            - cfg.phase_accumulation_reference_amplitude
        )
        return float(
            1.0
            + self._smoothstep(amplitude_position)
            * (cfg.phase_accumulation_high_amplitude_gain_ratio - 1.0)
        )

    def phase_accumulation_target_hardening_state(
        self, state: RIVASandState
    ) -> float:
        """Target amplitude-weighted activity for the next host increment."""
        return float(
            self.phase_accumulation_amplitude_gain(state)
            * self.phase_accumulation_target_activity(state)
        )

    def phase_accumulation_activity(self, state: RIVASandState) -> float:
        """Committed activity held constant through all constitutive substeps."""
        return float(
            max(getattr(state, "phase_accumulation_hardening_state", 0.0), 0.0)
        )

    def phase_accumulation_hardening_multiplier(
        self, state: RIVASandState
    ) -> float:
        return float(
            1.0
            + self.parameters.phase_accumulation_hardening_gain
            * self.phase_accumulation_activity(state)
        )

    def hardening_prefactor_for_state(
        self, pressure: float, state: RIVASandState
    ) -> float:
        hardening = super().hardening_prefactor_for_state(pressure, state)
        return float(
            hardening * self.phase_accumulation_hardening_multiplier(state)
        )

    def advance_fixed(
        self,
        initial: RIVASandAccumulationControlState,
        deps: Tensor,
        nsub: int = 1,
    ):
        """Advance backbone substeps with one committed shakedown update.

        The hardening state used by the local integrator is the value committed
        at the start of the host increment.  Its new target is accepted only
        after all fixed constitutive substeps finish.  Consequently, changing
        ``nsub`` refines the stress integration without changing the number of
        shakedown-memory updates.
        """
        state, info = super().advance_fixed(initial, deps, nsub)
        if not isinstance(state, RIVASandAccumulationControlState):
            raise TypeError("accumulation-control update lost its state type")
        state.phase_accumulation_hardening_state = (
            self.phase_accumulation_target_hardening_state(state)
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
            phase_accumulation_density_gate=self.phase_accumulation_density_gate(
                current
            ),
            phase_accumulation_bias_gate=self.phase_accumulation_bias_gate(current),
            phase_accumulation_memory=self.phase_accumulation_memory(current),
            phase_accumulation_onset=self.phase_accumulation_onset(current),
            phase_accumulation_activity=self.phase_accumulation_activity(current),
            phase_accumulation_target_activity=(
                self.phase_accumulation_target_activity(current)
            ),
            phase_accumulation_amplitude_gain=(
                self.phase_accumulation_amplitude_gain(current)
            ),
            phase_accumulation_target_hardening_state=(
                self.phase_accumulation_target_hardening_state(current)
            ),
            phase_accumulation_hardening_multiplier=(
                self.phase_accumulation_hardening_multiplier(current)
            ),
        )
        return values


__all__ = [
    "RIVASandAccumulationControlModel",
    "RIVASandAccumulationControlParameters",
    "RIVASandAccumulationControlState",
]
