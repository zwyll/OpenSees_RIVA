"""Loose-biased shear-flow successor to the signed-PT checkpoint.

This private prototype changes only the compact low-density, static-bias
window already used by :mod:`RIVASandLooseUnbiasedCorrection`.  It separates
the early-cycle shear response from a smoothly delayed plastic-work hardening
state while retaining the signed phase-volume law that removed the effective-
stress-path hourglass.

Setting ``loose_shear_flow_enabled=False`` recovers the parent checkpoint
exactly.  This is research code, not production RIVA-Sand.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from rivasand_port.model import RIVASandState, Tensor, invariants, tensor_norm

from RIVASandLooseUnbiasedCorrection import (
    RIVASandLooseUnbiasedCorrectionModel,
    RIVASandLooseUnbiasedCorrectionParameters,
    RIVASandLooseUnbiasedCorrectionState,
)


@dataclass(frozen=True)
class RIVASandLooseBiasedShearFlowParameters(
    RIVASandLooseUnbiasedCorrectionParameters
):
    """Signed-PT checkpoint plus isolated loose-biased shear evolution."""

    loose_phase_wave_activity_scale: float = 0.45
    loose_phase_mean_activity_scale: float = 1.40
    loose_shear_flow_enabled: bool = True
    loose_shear_cyclic_modulus_reduction: float = 0.85
    loose_shear_modulus_scale: float = 0.75
    loose_shear_early_hardening_multiplier: float = 1.40
    loose_shear_late_hardening_multiplier: float = 3.50
    loose_shear_hardening_memory_onset: float = 0.035
    loose_shear_hardening_memory_width: float = 0.010
    loose_shear_branch_compliance_gain: float = 1.0
    loose_shear_branch_compliance_onset: float = 0.60
    loose_shear_branch_compliance_full: float = 0.95
    loose_shear_ratchet_rate: float = 4.0
    loose_shear_ratchet_capacity: float = 0.0172
    loose_shear_ratchet_pressure_exponent: float = 1.0

    def __post_init__(self) -> None:
        super().__post_init__()
        if not 0.0 <= self.loose_shear_cyclic_modulus_reduction < 1.0:
            raise ValueError("loose shear modulus reduction must lie in [0,1)")
        if not 0.0 < self.loose_shear_modulus_scale <= 1.0:
            raise ValueError("loose shear modulus scale must lie in (0,1]")
        if self.loose_shear_early_hardening_multiplier <= 0.0:
            raise ValueError("loose early hardening multiplier must be positive")
        if (
            self.loose_shear_late_hardening_multiplier
            < self.loose_shear_early_hardening_multiplier
        ):
            raise ValueError("loose late hardening cannot be below early hardening")
        if self.loose_shear_hardening_memory_onset < 0.0:
            raise ValueError("loose hardening onset must be nonnegative")
        if self.loose_shear_hardening_memory_width <= 0.0:
            raise ValueError("loose hardening width must be positive")
        if self.loose_shear_branch_compliance_gain < 0.0:
            raise ValueError("loose branch compliance gain must be nonnegative")
        if not (
            0.0 <= self.loose_shear_branch_compliance_onset
            < self.loose_shear_branch_compliance_full <= 1.0
        ):
            raise ValueError("invalid loose branch compliance interval")
        if (
            self.loose_shear_ratchet_rate < 0.0
            or self.loose_shear_ratchet_capacity < 0.0
            or self.loose_shear_ratchet_pressure_exponent < 0.0
        ):
            raise ValueError("invalid loose ratchet controls")


@dataclass
class RIVASandLooseBiasedShearFlowState(
    RIVASandLooseUnbiasedCorrectionState
):
    """Adds a host-committed loose-sand plastic-work hardening state."""

    loose_shear_lambda_anchor: float = 0.0
    loose_shear_hardening_state: float = 0.0
    loose_shear_gate_value: float = 0.0

    def copy(self) -> "RIVASandLooseBiasedShearFlowState":
        values = {}
        for item in fields(RIVASandLooseBiasedShearFlowState):
            value = getattr(self, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        return RIVASandLooseBiasedShearFlowState(**values)


class RIVASandLooseBiasedShearFlowModel(
    RIVASandLooseUnbiasedCorrectionModel
):
    """Loose-biased loop correction with delayed committed hardening."""

    parameters: RIVASandLooseBiasedShearFlowParameters
    state: RIVASandLooseBiasedShearFlowState | None

    def __init__(
        self, parameters: RIVASandLooseBiasedShearFlowParameters | None = None
    ):
        super().__init__(parameters or RIVASandLooseBiasedShearFlowParameters())

    def initialize(
        self, stress: Tensor, void_ratio: float
    ) -> RIVASandLooseBiasedShearFlowState:
        base = super().initialize(stress, void_ratio)
        values = {
            item.name: getattr(base, item.name)
            for item in fields(RIVASandLooseUnbiasedCorrectionState)
        }
        for name, value in tuple(values.items()):
            if isinstance(value, np.ndarray):
                values[name] = value.copy()
        self.state = RIVASandLooseBiasedShearFlowState(**values)
        return self.state.copy()

    def begin_cyclic_phase(
        self, *, reference_stress: Tensor | None = None
    ) -> RIVASandLooseBiasedShearFlowState:
        super().begin_cyclic_phase(reference_stress=reference_stress)
        if not isinstance(self.state, RIVASandLooseBiasedShearFlowState):
            raise TypeError("loose shear-flow update lost its state type")
        self.state.loose_shear_lambda_anchor = float(self.state.lambda_total)
        self.state.loose_shear_hardening_state = 0.0
        self.state.loose_shear_gate_value = float(
            super().loose_phase_gate(self.state)
        )
        return self.state.copy()

    def loose_phase_gate(self, state: RIVASandState) -> float:
        """Return the density/static-bias gate frozen at cyclic activation."""
        if (
            isinstance(state, RIVASandLooseBiasedShearFlowState)
            and state.cyclic_phase_active
        ):
            return float(state.loose_shear_gate_value)
        return super().loose_phase_gate(state)

    def loose_shear_memory(self, state: RIVASandState) -> float:
        anchor = float(getattr(state, "loose_shear_lambda_anchor", 0.0))
        return float(max(state.lambda_total - anchor, 0.0))

    def loose_shear_target_hardening_state(self, state: RIVASandState) -> float:
        cfg = self.parameters
        if not cfg.loose_shear_flow_enabled or not state.cyclic_phase_active:
            return 0.0
        transition = self._smoothstep(
            (
                self.loose_shear_memory(state)
                - cfg.loose_shear_hardening_memory_onset
            )
            / cfg.loose_shear_hardening_memory_width
        )
        return float(self.loose_phase_gate(state) * transition)

    def loose_shear_hardening_multiplier(self, state: RIVASandState) -> float:
        cfg = self.parameters
        activity = float(max(
            getattr(state, "loose_shear_hardening_state", 0.0), 0.0
        ))
        return float(
            cfg.loose_shear_early_hardening_multiplier
            + activity
            * (
                cfg.loose_shear_late_hardening_multiplier
                - cfg.loose_shear_early_hardening_multiplier
            )
        )

    def hardening_prefactor_for_state(
        self, pressure: float, state: RIVASandState
    ) -> float:
        hardening = super().hardening_prefactor_for_state(pressure, state)
        cfg = self.parameters
        if not cfg.loose_shear_flow_enabled:
            return hardening
        gate = self.loose_phase_gate(state)
        if gate <= 1.0e-14 or not state.cyclic_phase_active:
            return hardening
        parent_multiplier = 1.0 + gate * (
            cfg.loose_stabilization_hardening_multiplier - 1.0
        )
        evolved_multiplier = self.loose_shear_hardening_multiplier(state)
        target_multiplier = 1.0 + gate * (evolved_multiplier - 1.0)
        return float(
            hardening
            / max(parent_multiplier, 1.0e-14)
            * target_multiplier
        )

    def cyclic_flow_factors(
        self, pressure: float, state: RIVASandState
    ) -> tuple[float, float]:
        shear, hardening = super().cyclic_flow_factors(pressure, state)
        cfg = self.parameters
        if (
            not cfg.loose_shear_flow_enabled
            or not cfg.cyclic_flow_correction_enabled
            or state.amplitude_reversals < cfg.cyclic_flow_minimum_reversals
        ):
            return shear, hardening
        gate = self.loose_phase_gate(state)
        if gate <= 1.0e-14:
            return shear, hardening
        pressure_ratio = np.clip(
            pressure / max(state.pressure_anchor, cfg.p_min), 0.0, 1.0
        )
        activity = pressure_ratio**cfg.cyclic_flow_pressure_exponent
        loose_shear = 1.0 - cfg.loose_shear_cyclic_modulus_reduction * activity
        return float(shear + gate * (loose_shear - shear)), hardening

    def moduli_for_state(
        self, pressure: float, state: RIVASandState
    ) -> tuple[float, float]:
        shear, bulk = super().moduli_for_state(pressure, state)
        cfg = self.parameters
        if not cfg.loose_shear_flow_enabled or not state.cyclic_phase_active:
            return shear, bulk
        gate = self.loose_phase_gate(state)
        scale = 1.0 + gate * (cfg.loose_shear_modulus_scale - 1.0)
        transition = self._smoothstep(
            (
                self.branch_progress(state)
                - cfg.loose_shear_branch_compliance_onset
            )
            / (
                cfg.loose_shear_branch_compliance_full
                - cfg.loose_shear_branch_compliance_onset
            )
        )
        branch = 1.0 / (
            1.0 + gate * cfg.loose_shear_branch_compliance_gain * transition
        )
        multiplier = scale * branch
        return float(shear * multiplier), float(bulk)

    def bias_ratchet_activity(self, state: RIVASandState) -> float:
        parent = super().bias_ratchet_activity(state)
        cfg = self.parameters
        if not cfg.loose_shear_flow_enabled:
            return parent
        gate = self.loose_phase_gate(state)
        loose = self.phase_amplitude_activity(state)
        return float(parent + gate * (loose - parent))

    def bias_ratchet_capacity(self, state: RIVASandState) -> float:
        parent = super().bias_ratchet_capacity(state)
        cfg = self.parameters
        if not cfg.loose_shear_flow_enabled:
            return parent
        gate = self.loose_phase_gate(state)
        loose = (
            cfg.loose_shear_ratchet_capacity
            * self.bias_ratchet_activity(state)
        )
        return float(parent + gate * (loose - parent))

    def bias_ratchet_rate_factor(self, state: RIVASandState) -> float:
        parent = super().bias_ratchet_rate_factor(state)
        cfg = self.parameters
        if not cfg.loose_shear_flow_enabled:
            return parent
        gate = self.loose_phase_gate(state)
        capacity = self.bias_ratchet_capacity(state)
        if capacity <= 1.0e-14 or state.bias_ratchet_strain >= capacity:
            loose = 0.0
        else:
            _, pressure, _ = invariants(state.stress)
            pressure_ratio = np.clip(
                pressure / max(state.pressure_anchor, cfg.p_min), 0.0, 1.0
            ) ** cfg.loose_shear_ratchet_pressure_exponent
            saturation = max(
                1.0 - state.bias_ratchet_strain / capacity, 0.0
            )
            loose = (
                cfg.loose_shear_ratchet_rate
                * self.bias_ratchet_activity(state)
                * pressure_ratio
                * saturation
            )
        return float(parent + gate * (loose - parent))

    def _bias_ratchet_increment(
        self, old: RIVASandState, trial: RIVASandState
    ) -> tuple[float, Tensor]:
        """Evaluate the isolated loose ratchet once per local update.

        The inherited update queries activity and capacity repeatedly.  This
        direct specialization retains the same equations while avoiding those
        repeated density/bias and invariant evaluations in the Python oracle.
        """
        cfg = self.parameters
        gate = self.loose_phase_gate(old)
        if not cfg.loose_shear_flow_enabled or gate <= 1.0e-14:
            return super()._bias_ratchet_increment(old, trial)
        parent_activity = super().bias_ratchet_activity(old)
        activity = parent_activity + gate * (
            self.phase_amplitude_activity(old) - parent_activity
        )
        bias = self.projected_bias(old)
        if bias <= 1.0e-14:
            parent_capacity = 0.0
        else:
            parent_capacity = (
                cfg.bias_ratchet_limit
                * activity
                * (
                    cfg.bias_ratchet_reference_bias / bias
                ) ** cfg.bias_ratchet_bias_exponent
            )
        loose_capacity = cfg.loose_shear_ratchet_capacity * activity
        capacity = parent_capacity + gate * (
            loose_capacity - parent_capacity
        )
        delta_lambda = max(trial.lambda_total - old.lambda_total, 0.0)
        direction_norm = tensor_norm(old.cyclic_direction)
        projection = float(np.sum(
            old.static_bias_tensor * old.cyclic_direction
        ))
        if (
            capacity <= 1.0e-14
            or delta_lambda <= 0.0
            or direction_norm <= 1.0e-12
            or old.bias_ratchet_strain >= capacity
        ):
            return 0.0, np.zeros((3, 3), dtype=float)
        _, pressure, _ = invariants(old.stress)
        pressure_ratio = np.clip(
            pressure / max(old.pressure_anchor, cfg.p_min), 0.0, 1.0
        ) ** cfg.loose_shear_ratchet_pressure_exponent
        saturation = max(
            1.0 - old.bias_ratchet_strain / capacity, 0.0
        )
        loose_rate = (
            cfg.loose_shear_ratchet_rate
            * activity
            * pressure_ratio
            * saturation
        )
        stress_level = max(
            old.pressure_anchor / cfg.bias_reference_pressure - 1.0, 0.0
        ) ** cfg.bias_ratchet_pressure_exponent
        bias_factor = (
            cfg.bias_ratchet_reference_bias / max(bias, 1.0e-14)
        ) ** cfg.bias_ratchet_bias_exponent
        parent_rate = (
            cfg.bias_ratchet_rate
            * stress_level
            * bias_factor
            * saturation
        )
        rate = parent_rate + gate * (loose_rate - parent_rate)
        increment = min(
            rate * delta_lambda,
            max(capacity - old.bias_ratchet_strain, 0.0),
        )
        direction = (
            np.sign(projection) * old.cyclic_direction / direction_norm
        )
        return float(increment), direction

    def advance_fixed(
        self,
        initial: RIVASandLooseBiasedShearFlowState,
        deps: Tensor,
        nsub: int = 1,
    ):
        state, info = super().advance_fixed(initial, deps, nsub)
        if not isinstance(state, RIVASandLooseBiasedShearFlowState):
            raise TypeError("loose shear-flow update lost its state type")
        state.loose_shear_hardening_state = (
            self.loose_shear_target_hardening_state(state)
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
            loose_shear_memory=self.loose_shear_memory(current),
            loose_shear_hardening_state=float(getattr(
                current, "loose_shear_hardening_state", 0.0
            )),
            loose_shear_gate_value=float(getattr(
                current, "loose_shear_gate_value", 0.0
            )),
            loose_shear_hardening_multiplier=(
                self.loose_shear_hardening_multiplier(current)
            ),
            loose_shear_branch_compliance_multiplier=(
                1.0
                / (
                    1.0
                    + self.loose_phase_gate(current)
                    * self.parameters.loose_shear_branch_compliance_gain
                    * self._smoothstep(
                        (
                            self.branch_progress(current)
                            - self.parameters.loose_shear_branch_compliance_onset
                        )
                        / (
                            self.parameters.loose_shear_branch_compliance_full
                            - self.parameters.loose_shear_branch_compliance_onset
                        )
                    )
                )
            ),
        )
        return values


__all__ = [
    "RIVASandLooseBiasedShearFlowModel",
    "RIVASandLooseBiasedShearFlowParameters",
    "RIVASandLooseBiasedShearFlowState",
]
