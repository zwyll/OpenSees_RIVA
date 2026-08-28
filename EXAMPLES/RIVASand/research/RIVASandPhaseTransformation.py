"""Restrained phase-transformation research successor to production RIVA-Sand.

This private prototype starts directly from the frozen RIVA-Sand oracle.
One smooth transformation coordinate controls both
the within-branch shear compliance and compatible volumetric flow. The
volumetric response is integrated incrementally; it never assigns pressure or
pore-pressure ratio directly.

Setting ``phase_transformation_enabled=False`` recovers production RIVA-Sand
exactly. This is research code, not a production release.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from rivasand_port.model import (
    I3,
    IntegrationInfo,
    RIVASandModel,
    RIVASandParameters,
    RIVASandState,
    SQRT23,
    Tensor,
    invariants,
    tensor_norm,
)


@dataclass(frozen=True)
class RIVASandPhaseTransformationParameters(RIVASandParameters):
    """Frozen RIVA-Sand plus one coupled transformation-zone law."""

    phase_transformation_enabled: bool = True

    # Dense biased shear response. These controls are inactive outside the
    # phase branch, so disabling the branch recovers production RIVA-Sand.
    phase_cyclic_shear_modulus_reduction: float = 0.813
    branch_compliance_enabled: bool = True
    branch_compliance_bias_reference: float = 0.5384061785684372
    branch_compliance_bias_exponent: float = 1.0
    branch_compliance_minimum: float = 0.10

    # Smooth stiff-to-soft transition within each half-cycle. The transition
    # remains monotonic so reversal resets retain finite hysteretic area.
    phase_compliance_peak: float = 6.0
    phase_compliance_bell_gain: float = 3.5
    phase_memory_hardening_reduction: float = 2.8
    phase_memory_reference_volume: float = 5.0e-5
    phase_bias_hardening_intercept: float = 87.291302
    phase_bias_hardening_exponent: float = 3.2
    phase_compliance_shape: float = 1.55
    phase_compliance_location: float = 0.442
    phase_compliance_half_width: float = 0.787
    branch_directional_balance: float = 0.050
    branch_balance_bias_exponent: float = 4.0
    branch_balance_bias_cap: float = 3.9

    # Phase-transformation stress ratio and smoothness. The irreversible rate
    # multiplies deviatoric strain path length in the contractile region. The
    # reversible scale is a bounded volumetric-strain potential driven by one
    # fixed signed static-bias direction and approached over a strain length.
    phase_ratio: float = 0.62
    phase_width: float = 0.50
    phase_contraction_rate: float = 4.0e-4
    phase_reversible_scale: float = -3.00e-3
    phase_reversible_mean_scale: float = 6.50e-4
    phase_reversible_relaxation_strain: float = 0.0048
    phase_reversible_relaxation_bias_exponent: float = 1.6
    phase_potential_anchor_fraction: float = 1.0
    phase_pressure_exponent: float = 1.6
    phase_wave_pressure_exponent: float = 1.0
    phase_mean_pressure_exponent: float = 0.0
    phase_amplitude_onset_ratio: float = 0.55
    phase_amplitude_full_ratio: float = 0.90
    phase_volume_density_onset: float = 0.53
    phase_volume_density_full: float = 0.90
    phase_volume_replacement_density_full: float = 0.66
    phase_intermediate_wave_multiplier: float = 1.1333333333333333
    phase_intermediate_mean_multiplier: float = 10.76923076923077
    phase_intermediate_relaxation_ratio: float = 0.10416666666666667
    phase_bias_exponent: float = 2.0
    phase_wave_bias_exponent: float = -0.25
    phase_mean_bias_exponent: float = 12.0

    def __post_init__(self) -> None:
        super().__post_init__()
        if not 0.0 <= self.phase_cyclic_shear_modulus_reduction < 1.0:
            raise ValueError("phase cyclic shear reduction must lie in [0,1)")
        if self.branch_compliance_bias_reference <= 0.0:
            raise ValueError("branch compliance bias reference must be positive")
        if self.branch_compliance_bias_exponent < 0.0:
            raise ValueError("branch compliance bias exponent must be nonnegative")
        if not 0.0 < self.branch_compliance_minimum <= 1.0:
            raise ValueError("branch compliance minimum must lie in (0,1]")
        if self.phase_compliance_peak < 0.0:
            raise ValueError("phase compliance peak must be nonnegative")
        if self.phase_compliance_bell_gain < 0.0:
            raise ValueError("phase compliance bell gain must be nonnegative")
        if self.phase_memory_hardening_reduction < 0.0:
            raise ValueError("phase memory hardening reduction must be nonnegative")
        if self.phase_memory_reference_volume <= 0.0:
            raise ValueError("phase memory reference volume must be positive")
        if self.phase_bias_hardening_intercept < 0.0:
            raise ValueError("phase bias hardening intercept must be nonnegative")
        if self.phase_bias_hardening_exponent <= 0.0:
            raise ValueError("phase bias hardening exponent must be positive")
        if self.branch_balance_bias_cap <= 0.0:
            raise ValueError("branch balance bias cap must be positive")
        if self.phase_compliance_shape <= 0.0:
            raise ValueError("phase compliance shape must be positive")
        if not 0.0 < self.phase_compliance_location < 1.0:
            raise ValueError("phase compliance location must lie in (0,1)")
        if not 0.0 < self.phase_compliance_half_width <= 1.0:
            raise ValueError("phase compliance half width must lie in (0,1]")
        if self.phase_ratio <= 0.0 or self.phase_width <= 0.0:
            raise ValueError("phase ratio and width must be positive")
        if self.phase_contraction_rate < 0.0:
            raise ValueError("phase contraction rate must be nonnegative")
        if not np.isfinite(self.phase_reversible_scale):
            raise ValueError("phase reversible scale must be finite")
        if not np.isfinite(self.phase_reversible_mean_scale):
            raise ValueError("phase reversible mean scale must be finite")
        if self.phase_reversible_relaxation_strain <= 0.0:
            raise ValueError("phase reversible relaxation strain must be positive")
        if self.phase_reversible_relaxation_bias_exponent < 0.0:
            raise ValueError("phase relaxation bias exponent must be nonnegative")
        if not 0.0 <= self.phase_potential_anchor_fraction <= 2.0:
            raise ValueError("phase potential anchor fraction must lie in [0,2]")
        if self.phase_pressure_exponent < 0.0:
            raise ValueError("phase pressure exponent must be nonnegative")
        if self.phase_wave_pressure_exponent < 0.0:
            raise ValueError("phase wave pressure exponent must be nonnegative")
        if self.phase_mean_pressure_exponent < 0.0:
            raise ValueError("phase mean pressure exponent must be nonnegative")
        if not (
            0.0 <= self.phase_amplitude_onset_ratio
            < self.phase_amplitude_full_ratio
        ):
            raise ValueError("invalid phase amplitude activation interval")
        if not (
            0.0 <= self.phase_volume_density_onset
            < self.phase_volume_density_full <= 1.0
        ):
            raise ValueError("invalid phase volume density interval")
        if not (
            self.phase_volume_density_onset
            < self.phase_volume_replacement_density_full <= 1.0
        ):
            raise ValueError("invalid phase volume replacement interval")
        if (
            self.phase_intermediate_wave_multiplier <= 0.0
            or self.phase_intermediate_mean_multiplier <= 0.0
            or self.phase_intermediate_relaxation_ratio <= 0.0
        ):
            raise ValueError("intermediate-density PT multipliers must be positive")
        if self.phase_bias_exponent < 0.0:
            raise ValueError("phase bias exponent must be nonnegative")
        if not np.isfinite(self.phase_wave_bias_exponent):
            raise ValueError("phase wave bias exponent must be finite")
        if self.phase_mean_bias_exponent < 0.0:
            raise ValueError("phase mean bias exponent must be nonnegative")


@dataclass
class RIVASandPhaseTransformationState(RIVASandState):
    """Three scalar memories required by the continuous volume evolution."""

    phase_irreversible_volume: float = 0.0
    phase_reversible_volume: float = 0.0
    phase_potential_anchor: float = 0.0

    def copy(self) -> "RIVASandPhaseTransformationState":
        values = {}
        for item in fields(RIVASandPhaseTransformationState):
            value = getattr(self, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        return RIVASandPhaseTransformationState(**values)


class RIVASandPhaseTransformationModel(RIVASandModel):
    """RIVA-Sand research model with coupled shear/volume transformation."""

    parameters: RIVASandPhaseTransformationParameters
    state: RIVASandPhaseTransformationState | None

    def __init__(
        self, parameters: RIVASandPhaseTransformationParameters | None = None
    ):
        super().__init__(parameters or RIVASandPhaseTransformationParameters())
        self._inside_host_substeps = False

    def initialize(
        self, stress: Tensor, void_ratio: float
    ) -> RIVASandPhaseTransformationState:
        base = super().initialize(stress, void_ratio)
        values = {
            item.name: getattr(base, item.name)
            for item in fields(RIVASandState)
        }
        for name, value in tuple(values.items()):
            if isinstance(value, np.ndarray):
                values[name] = value.copy()
        self.state = RIVASandPhaseTransformationState(**values)
        return self.state.copy()

    def begin_cyclic_phase(
        self, *, reference_stress: Tensor | None = None
    ) -> RIVASandPhaseTransformationState:
        super().begin_cyclic_phase(reference_stress=reference_stress)
        if not isinstance(self.state, RIVASandPhaseTransformationState):
            raise TypeError("phase-transformation update lost its state type")
        self.state.phase_irreversible_volume = 0.0
        self.state.phase_reversible_volume = 0.0
        # Starting volume remains zero and the relaxation law admits the first
        # dynamic increment continuously.  A bounded fraction of the static
        # potential permits calibration of mean pressure independently from
        # the cyclic pressure-wave amplitude.
        self.state.phase_potential_anchor = float(
            self.parameters.phase_potential_anchor_fraction
            * self.signed_phase_potential(self.state)
        )
        return self.state.copy()

    @staticmethod
    def _smoothstep(value: float) -> float:
        # This helper is called only with scalar coordinates, often dozens of
        # times in one constitutive substep.  ``np.clip`` dispatches through
        # NumPy's array protocol even for a scalar; direct bounds preserve the
        # same map without that repeated machinery.
        value = float(value)
        if value <= 0.0:
            return 0.0
        if value >= 1.0:
            return 1.0
        return value * value * (3.0 - 2.0 * value)

    def branch_progress(self, state: RIVASandState) -> float:
        """Normalized deviatoric excursion from the current reversal point."""
        cfg = self.parameters
        if (
            not cfg.phase_transformation_enabled
            or not cfg.branch_compliance_enabled
            or not state.cyclic_phase_active
            or state.amplitude_reversals < 1
            or state.cyclic_amplitude <= 1.0e-14
        ):
            return 0.0
        deviator, _, _ = invariants(state.stress)
        excursion = tensor_norm(deviator - state.last_reversal_deviator)
        span = 2.0 * max(
            state.pressure_anchor, cfg.p_min
        ) * state.cyclic_amplitude
        return float(np.clip(
            excursion / max(span, 1.0e-14), 0.0, 1.0
        ))

    def transformation_progress(self, state: RIVASandState) -> float:
        """Monotonic smooth transition from the stiff to soft branch state."""
        cfg = self.parameters
        if not cfg.phase_transformation_enabled:
            return 0.0
        progress = self.branch_progress(state)
        lower = cfg.phase_compliance_location - cfg.phase_compliance_half_width
        upper = cfg.phase_compliance_location + cfg.phase_compliance_half_width
        if progress <= lower:
            return 0.0
        if progress >= upper:
            return 1.0
        coordinate = (progress - lower) / max(upper - lower, 1.0e-12)
        return float(self._smoothstep(coordinate) ** cfg.phase_compliance_shape)

    def transformation_zone(self, state: RIVASandState) -> float:
        """Fixed-direction phase zone centered on low signed shear stress."""
        cfg = self.parameters
        if not cfg.phase_transformation_enabled:
            return 0.0
        potential = self.signed_phase_potential(state)
        return float(max(1.0 - potential * potential, 0.0))

    def branch_stress_phase(self, state: RIVASandState) -> float:
        """Signed dynamic stress coordinate about the admitted static state."""
        direction_norm = tensor_norm(state.cyclic_direction)
        if (
            direction_norm <= 1.0e-14
            or state.cyclic_amplitude <= 1.0e-14
        ):
            return 1.0
        direction = state.cyclic_direction / direction_norm
        deviator, _, _ = invariants(state.stress)
        static_deviator = state.geostatic_deviator + (
            state.pressure_anchor * state.static_bias_tensor
        )
        projection = float(np.sum((deviator - static_deviator) * direction))
        scale = state.pressure_anchor * state.cyclic_amplitude
        return float(np.clip(projection / max(scale, 1.0e-14), -1.0, 1.0))

    def phase_amplitude_activity(self, state: RIVASandState) -> float:
        cfg = self.parameters
        ratio = state.cyclic_amplitude / cfg.cyclic_amplitude_reference
        position = (
            ratio - cfg.phase_amplitude_onset_ratio
        ) / (
            cfg.phase_amplitude_full_ratio - cfg.phase_amplitude_onset_ratio
        )
        return self._smoothstep(position)

    def compliance_bias_activity(self, state: RIVASandState) -> float:
        cfg = self.parameters
        bias = self.projected_bias(state)
        if bias <= 1.0e-14:
            return 0.0
        return float(min(
            (bias / cfg.branch_compliance_bias_reference)
            ** cfg.branch_compliance_bias_exponent,
            2.0,
        ))

    def phase_flow_bias_activity(
        self, state: RIVASandState, *, exponent: float | None = None
    ) -> float:
        """Smooth inverse-bias scaling for transformation volume."""
        cfg = self.parameters
        bias = self.projected_bias(state)
        reference = cfg.branch_compliance_bias_reference
        if bias <= 1.0e-14:
            return 0.0
        onset = 0.25 * reference
        gate = self._smoothstep(bias / onset)
        power = cfg.phase_bias_exponent if exponent is None else exponent
        return float(min(
            gate * (reference / bias) ** power,
            2.0,
        ))

    def phase_dense_bias_gate(self, state: RIVASandState) -> float:
        """Activation shared by the dense biased shear and hardening laws."""
        bias = self.projected_bias(state)
        onset = 0.50 * self.parameters.branch_compliance_bias_reference
        return float(self.dense_state_weight(state) * self._smoothstep(
            bias / max(onset, 1.0e-14)
        ))

    def phase_volume_density_weight(self, state: RIVASandState) -> float:
        """Smoothly extend compatible PT volume below the dense shear gate."""
        cfg = self.parameters
        position = (
            self.initial_relative_density(state) - cfg.phase_volume_density_onset
        ) / (
            cfg.phase_volume_density_full - cfg.phase_volume_density_onset
        )
        return self._smoothstep(position)

    def phase_volume_gate(self, state: RIVASandState) -> float:
        """Density--bias amplitude of the fixed-direction PT volume law."""
        bias = self.projected_bias(state)
        onset = 0.50 * self.parameters.branch_compliance_bias_reference
        return float(self.phase_volume_density_weight(state) * self._smoothstep(
            bias / max(onset, 1.0e-14)
        ))

    def phase_volume_replacement_gate(self, state: RIVASandState) -> float:
        """Remove the branch-relative wave before its hourglass fold appears."""
        cfg = self.parameters
        density_position = (
            self.initial_relative_density(state) - cfg.phase_volume_density_onset
        ) / (
            cfg.phase_volume_replacement_density_full
            - cfg.phase_volume_density_onset
        )
        bias = self.projected_bias(state)
        bias_onset = 0.50 * cfg.branch_compliance_bias_reference
        return float(
            self._smoothstep(density_position)
            * self._smoothstep(bias / max(bias_onset, 1.0e-14))
        )

    def phase_volume_dense_transition(self, state: RIVASandState) -> float:
        """Blend intermediate PT controls into the frozen dense calibration."""
        cfg = self.parameters
        reference = self.reference_relative_density()
        position = (
            self.initial_relative_density(state) - reference
        ) / max(cfg.phase_volume_density_full - reference, 1.0e-14)
        return self._smoothstep(position)

    def branch_compliance_multiplier(self, state: RIVASandState) -> float:
        """Localize the checkpoint's reversal-induced compliance transition."""
        cfg = self.parameters
        if not cfg.phase_transformation_enabled:
            return 1.0
        density = self.dense_state_weight(state)
        bias = self.compliance_bias_activity(state)
        transition = self.transformation_progress(state)
        if density <= 1.0e-14 or bias <= 1.0e-14 or transition <= 1.0e-14:
            return 1.0
        direction_norm = tensor_norm(state.cyclic_direction)
        projection = float(
            np.sum(state.static_bias_tensor * state.cyclic_direction)
            / max(direction_norm, 1.0e-14)
        )
        signed_bias = float(
            np.sign(projection)
            * min(
                (
                    abs(projection) / cfg.branch_compliance_bias_reference
                ) ** cfg.branch_balance_bias_exponent,
                cfg.branch_balance_bias_cap,
            )
        )
        balance = np.exp(cfg.branch_directional_balance * signed_bias)
        compliance = density * bias * (
            cfg.phase_compliance_peak * transition * balance
            + cfg.phase_compliance_bell_gain * self.transformation_zone(state)
        )
        return float(np.clip(
            1.0 / (1.0 + compliance),
            cfg.branch_compliance_minimum,
            1.0,
        ))

    def cyclic_flow_factors(
        self, pressure: float, state: RIVASandState
    ) -> tuple[float, float]:
        """Blend the calibrated dense PT cut into the frozen flow factors."""
        shear, hardening = super().cyclic_flow_factors(pressure, state)
        cfg = self.parameters
        if (
            not cfg.phase_transformation_enabled
            or not cfg.cyclic_flow_correction_enabled
            or state.amplitude_reversals < cfg.cyclic_flow_minimum_reversals
        ):
            return shear, hardening
        pressure_ratio = np.clip(
            pressure / max(state.pressure_anchor, cfg.p_min), 0.0, 1.0
        )
        activity = pressure_ratio**cfg.cyclic_flow_pressure_exponent
        phase_shear = 1.0 - cfg.phase_cyclic_shear_modulus_reduction * activity
        gate = self.phase_dense_bias_gate(state)
        return float(shear + gate * (phase_shear - shear)), hardening

    def moduli_for_state(
        self, pressure: float, state: RIVASandState
    ) -> tuple[float, float]:
        shear, bulk = super().moduli_for_state(pressure, state)
        shear *= self.branch_compliance_multiplier(state)
        return float(shear), float(bulk)

    def bias_hardening_boost(self, state: RIVASandState) -> float:
        """Evolve excess static-bias hardening inside the PT zone.

        The frozen model's excess static-bias hardening is retained outside
        the transformation zone.  The optional memory term degrades only the
        dense, high-bias branch as compatible irreversible phase volume
        accumulates.  Its smooth buildup permits delayed cyclic mobility
        without an instantaneous cycle-count or amplitude knee.
        """
        frozen_boost = super().bias_hardening_boost(state)
        cfg = self.parameters
        if not cfg.phase_transformation_enabled:
            return frozen_boost
        if (
            not cfg.static_bias_enabled
            or not state.cyclic_phase_active
            or state.amplitude_reversals < cfg.bias_minimum_reversals
        ):
            return frozen_boost
        bias = self.projected_bias(state)
        amplitude = cfg.bias_amplitude_ratio * state.cyclic_amplitude
        margin = max(bias - amplitude, 0.0)
        crossing = max(amplitude / max(bias, 1.0e-14) - 1.0, 0.0)
        phase_boost = (
            cfg.phase_bias_hardening_intercept
            * bias**cfg.phase_bias_hardening_exponent
            * np.exp(-cfg.bias_crossing_decay * crossing)
            + cfg.bias_hardening_scale * margin**cfg.bias_margin_exponent
        )
        gate = self.phase_dense_bias_gate(state)
        boost = frozen_boost + gate * (phase_boost - frozen_boost)
        bias_position = (
            bias / cfg.branch_compliance_bias_reference - 1.0
        ) / 0.50
        high_bias_gate = self._smoothstep(bias_position)
        accumulated_volume = max(
            -float(getattr(state, "phase_irreversible_volume", 0.0)), 0.0
        )
        memory = 1.0 - np.exp(
            -accumulated_volume / cfg.phase_memory_reference_volume
        )
        log_multiplier = (
            -cfg.phase_memory_hardening_reduction
            * self.dense_state_weight(state)
            * self.transformation_zone(state)
            * high_bias_gate
            * memory
        )
        return float(boost * np.exp(log_multiplier))

    def phase_coordinates(
        self, state: RIVASandState
    ) -> tuple[float, float, float]:
        """Return loading stress ratio, PT ratio, and smooth phase sign."""
        cfg = self.parameters
        direction_norm = tensor_norm(state.cyclic_direction)
        if direction_norm <= 1.0e-14:
            return 0.0, 0.0, -1.0
        loading_direction = -state.cyclic_direction / direction_norm
        loading_eta = float(np.sum(state.alpha * loading_direction))
        _, pressure, _ = invariants(state.stress)
        _, md, _ = self.surfaces(pressure, state.void_ratio)
        eta_pt = cfg.phase_ratio * SQRT23 * md
        signed_phase = float(np.tanh(
            (loading_eta - eta_pt) / max(cfg.phase_width * eta_pt, 1.0e-12)
        ))
        return loading_eta, float(eta_pt), signed_phase

    def signed_phase_potential(self, state: RIVASandState) -> float:
        """Smooth potential measured along one fixed static-bias direction."""
        direction_norm = tensor_norm(state.static_bias_tensor)
        if direction_norm <= 1.0e-14:
            return 0.0
        direction = state.static_bias_tensor / direction_norm
        signed_eta = float(np.sum(state.alpha * direction))
        _, pressure, _ = invariants(state.stress)
        _, md, _ = self.surfaces(pressure, state.void_ratio)
        eta_pt = self.parameters.phase_ratio * SQRT23 * md
        return float(np.tanh(
            signed_eta / max(self.parameters.phase_width * eta_pt, 1.0e-12)
        ))

    def phase_activity(
        self,
        state: RIVASandState,
        *,
        bias_exponent: float | None = None,
        pressure_exponent: float | None = None,
    ) -> float:
        if not state.cyclic_phase_active or state.amplitude_reversals < 1:
            return 0.0
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
        return float(
            self.phase_volume_density_weight(state)
            * self.phase_flow_bias_activity(state, exponent=bias_exponent)
            * self.phase_amplitude_activity(state)
            * pressure_ratio
        )

    def phase_reversible_relaxation_multiplier(
        self, state: RIVASandState
    ) -> float:
        """Hook for research successors; unity preserves this checkpoint."""
        return 1.0

    def bias_reversible_volume_target(self, state: RIVASandState) -> float:
        """Retain the frozen wave outside the dense transformation branch."""
        target = super().bias_reversible_volume_target(state)
        if not self.parameters.phase_transformation_enabled:
            return target
        target *= 1.0 - self.phase_volume_replacement_gate(state)
        # The frozen update removes the prior reversible bias before each
        # constitutive substep. Reapply the committed PT volume as its target
        # so no additional pressure rebuild is needed inside the substep loop.
        if self._inside_host_substeps:
            target += float(getattr(state, "phase_reversible_volume", 0.0))
        return float(target)

    def _rebuild_phase_compatible_pressure(
        self, state: RIVASandPhaseTransformationState
    ) -> None:
        effective_total = state.physical_eps_v_total + state.bias_reversible_volume
        state.eps_v_total = float(effective_total)
        state.eps_v_confining = float(
            state.eps_v_total
            - state.eps_v_irreversible
            - state.eps_v_reversible
        )
        deviator, _, _ = invariants(state.stress)
        pressure, at_floor = self.pressure_from_confining_strain(
            state.eps_v_confining, state.pressure_anchor
        )
        state.alpha = deviator / pressure
        state.stress = deviator - pressure * I3
        state.pressure_floor_hits += int(at_floor)

    def _apply_host_phase_volume(
        self,
        old: RIVASandPhaseTransformationState,
        state: RIVASandPhaseTransformationState,
        deps: Tensor,
    ) -> RIVASandPhaseTransformationState:
        """Advance phase memory exactly once for one physical host increment."""
        cfg = self.parameters
        deps_v = float(np.trace(deps))
        deps_d = np.asarray(deps, dtype=float) - deps_v * I3 / 3.0
        path_increment = tensor_norm(deps_d)
        activity = 0.5 * (self.phase_activity(old) + self.phase_activity(state))
        old_wave_activity = self.phase_activity(
            old,
            bias_exponent=cfg.phase_wave_bias_exponent,
            pressure_exponent=cfg.phase_wave_pressure_exponent,
        )
        new_wave_activity = self.phase_activity(
            state,
            bias_exponent=cfg.phase_wave_bias_exponent,
            pressure_exponent=cfg.phase_wave_pressure_exponent,
        )
        old_mean_activity = self.phase_activity(
            old,
            bias_exponent=cfg.phase_mean_bias_exponent,
            pressure_exponent=cfg.phase_mean_pressure_exponent,
        )
        new_mean_activity = self.phase_activity(
            state,
            bias_exponent=cfg.phase_mean_bias_exponent,
            pressure_exponent=cfg.phase_mean_pressure_exponent,
        )
        _, _, signed_phase = self.phase_coordinates(state)
        contractile = 0.5 * (1.0 - signed_phase)
        irreversible_increment = (
            -cfg.phase_contraction_rate
            * activity
            * self.transformation_zone(state)
            * contractile
            * path_increment
        )
        state.phase_irreversible_volume = float(
            old.phase_irreversible_volume + irreversible_increment
        )
        state.eps_v_irreversible = float(
            state.eps_v_irreversible + irreversible_increment
        )
        old_potential = self.signed_phase_potential(old)
        new_potential = self.signed_phase_potential(state)
        old_phase_memory = 1.0 - np.exp(
            -max(-old.phase_irreversible_volume, 0.0)
            / cfg.phase_memory_reference_volume
        )
        new_phase_memory = 1.0 - np.exp(
            -max(-state.phase_irreversible_volume, 0.0)
            / cfg.phase_memory_reference_volume
        )
        dense_transition = 0.5 * (
            self.phase_volume_dense_transition(old)
            + self.phase_volume_dense_transition(state)
        )
        wave_multiplier = (
            cfg.phase_intermediate_wave_multiplier
            + dense_transition * (1.0 - cfg.phase_intermediate_wave_multiplier)
        )
        mean_multiplier = (
            cfg.phase_intermediate_mean_multiplier
            + dense_transition * (1.0 - cfg.phase_intermediate_mean_multiplier)
        )
        old_reversible_target = (
            cfg.phase_reversible_scale
            * wave_multiplier
            * old_wave_activity
            * (old_potential - old.phase_potential_anchor)
            + cfg.phase_reversible_mean_scale
            * mean_multiplier
            * old_mean_activity
            * old_phase_memory
        )
        new_reversible_target = (
            cfg.phase_reversible_scale
            * wave_multiplier
            * new_wave_activity
            * (new_potential - state.phase_potential_anchor)
            + cfg.phase_reversible_mean_scale
            * mean_multiplier
            * new_mean_activity
            * new_phase_memory
        )
        reversible_target = 0.5 * (
            old_reversible_target + new_reversible_target
        )
        bias_ratio = max(
            self.projected_bias(state) / cfg.branch_compliance_bias_reference,
            1.0,
        )
        relaxation_strain = (
            cfg.phase_reversible_relaxation_strain
            * (
                cfg.phase_intermediate_relaxation_ratio
                + dense_transition
                * (1.0 - cfg.phase_intermediate_relaxation_ratio)
            )
            * bias_ratio**cfg.phase_reversible_relaxation_bias_exponent
        )
        relaxation_multiplier = max(
            self.phase_reversible_relaxation_multiplier(state), 0.0
        )
        relaxation = 1.0 - np.exp(
            -relaxation_multiplier * path_increment / relaxation_strain
        )
        state.phase_reversible_volume = float(
            old.phase_reversible_volume
            + relaxation * (reversible_target - old.phase_reversible_volume)
        )
        # The internal substep loop already restored the committed phase
        # volume. Replace it, rather than adding the new value on top.
        state.bias_reversible_volume = float(
            state.bias_reversible_volume
            - old.phase_reversible_volume
            + state.phase_reversible_volume
        )
        self._rebuild_phase_compatible_pressure(state)
        return state

    def _forward_euler(
        self,
        old: RIVASandPhaseTransformationState,
        deps: Tensor,
        *,
        force_reversal: bool = False,
        allow_legacy_reversal: bool = True,
    ) -> RIVASandPhaseTransformationState:
        state = super()._forward_euler(
            old,
            deps,
            force_reversal=force_reversal,
            allow_legacy_reversal=allow_legacy_reversal,
        )
        if not isinstance(state, RIVASandPhaseTransformationState):
            raise TypeError("phase-transformation update lost its state type")
        return state

    def advance_fixed(
        self,
        initial: RIVASandPhaseTransformationState,
        deps: Tensor,
        nsub: int = 1,
    ) -> tuple[RIVASandPhaseTransformationState, IntegrationInfo]:
        """Refine the backbone while advancing phase memory once per host step."""
        if (
            not self.parameters.phase_transformation_enabled
            or (
                self.phase_dense_bias_gate(initial) <= 1.0e-14
                and self.phase_volume_gate(initial) <= 1.0e-14
            )
        ):
            return super().advance_fixed(initial, deps, nsub)
        self._inside_host_substeps = True
        try:
            state, info = super().advance_fixed(initial, deps, nsub)
        finally:
            self._inside_host_substeps = False
        return self._apply_host_phase_volume(initial, state, deps), info

    def dss_history_values(
        self, state: RIVASandState | None = None
    ) -> dict[str, float]:
        current = state or self.state
        if current is None:
            raise RuntimeError("initialize the model first")
        values = super().dss_history_values(current)
        loading_eta, eta_pt, signed_phase = self.phase_coordinates(current)
        values.update(
            lambda_total=current.lambda_total,
            ep_eq_since_reversal=current.ep_eq_since_reversal,
            branch_progress=self.branch_progress(current),
            branch_compliance_multiplier=self.branch_compliance_multiplier(current),
            branch_bias_sign=float(np.sign(
                np.sum(current.static_bias_tensor * current.cyclic_direction)
            )),
            phase_dense_bias_gate=self.phase_dense_bias_gate(current),
            phase_volume_density_weight=self.phase_volume_density_weight(current),
            phase_volume_gate=self.phase_volume_gate(current),
            phase_volume_replacement_gate=self.phase_volume_replacement_gate(current),
            phase_volume_dense_transition=self.phase_volume_dense_transition(current),
            transformation_progress=self.transformation_progress(current),
            transformation_zone=self.transformation_zone(current),
            branch_stress_phase=self.branch_stress_phase(current),
            phase_loading_eta=loading_eta,
            phase_transformation_eta=eta_pt,
            phase_signed_state=signed_phase,
            phase_activity=self.phase_activity(current),
            phase_wave_activity=self.phase_activity(
                current,
                bias_exponent=self.parameters.phase_wave_bias_exponent,
                pressure_exponent=self.parameters.phase_wave_pressure_exponent,
            ),
            phase_mean_activity=self.phase_activity(
                current,
                bias_exponent=self.parameters.phase_mean_bias_exponent,
                pressure_exponent=self.parameters.phase_mean_pressure_exponent,
            ),
            phase_potential=self.signed_phase_potential(current),
            phase_irreversible_volume=float(getattr(
                current, "phase_irreversible_volume", 0.0
            )),
            phase_reversible_volume=float(getattr(
                current, "phase_reversible_volume", 0.0
            )),
        )
        return values


__all__ = [
    "RIVASandPhaseTransformationModel",
    "RIVASandPhaseTransformationParameters",
    "RIVASandPhaseTransformationState",
]
