"""Directional mapping-backstress successor to the qualified DSS checkpoint.

This private research kernel replaces the fragile high-static-bias scalar
overlay with a tensorial mapping center and directional fabric.  The new
mechanism is deliberately isolated from the calibrated RIVA-Sand histories:
it is inactive at and below the alpha=0.15 loose-biased calibration and for
the intermediate/dense Ottawa tests.

The active update uses a fixed, stress-corrected midpoint iteration.  The
backstress rate contributes directly to the plastic denominator, and the
plastic multiplier is bounded so the incremental deviatoric stress cannot do
negative work in the imposed strain direction.  No adaptive substepping is
used.  This is a research successor, not production RIVA-Sand.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields

import numpy as np

from rivasand_port.model import (
    FAC23,
    I3,
    SQRT23,
    RIVASandState,
    Tensor,
    invariants,
    tensor_norm,
)

from RIVASandLooseBiasedShearFlow import (
    RIVASandLooseBiasedShearFlowModel,
    RIVASandLooseBiasedShearFlowParameters,
    RIVASandLooseBiasedShearFlowState,
)


def _unit(value: Tensor) -> Tensor:
    norm = tensor_norm(value)
    if norm <= 1.0e-14:
        return np.zeros((3, 3), dtype=float)
    return np.asarray(value, dtype=float) / norm


@dataclass(frozen=True)
class RIVASandMappingBackstressParameters(
    RIVASandLooseBiasedShearFlowParameters
):
    """Controls for the isolated directional mapping correction."""

    mapping_backstress_enabled: bool = True
    mapping_bias_onset: float = 0.34
    mapping_bias_full: float = 0.54
    mapping_density_full: float = 0.55
    mapping_density_cutoff: float = 0.64

    mapping_backstress_rate: float = 11000.0
    mapping_backstress_capacity_fraction: float = 0.90
    mapping_center_limit_ratio: float = 0.97
    mapping_core_radius_ratio: float = 0.020
    mapping_ray_flow_weight: float = 0.16

    mapping_fabric_dilation_rate: float = 120.0
    mapping_fabric_recovery_rate: float = 30.0
    mapping_fabric_saturation: float = 0.75
    mapping_fabric_flow_weight: float = 0.10
    mapping_fabric_dilatancy_weight: float = 0.20
    mapping_directional_ratchet_weight: float = 0.020
    mapping_high_bias_contraction_onset: float = 0.62
    mapping_high_bias_contraction_full: float = 0.82
    mapping_high_bias_contraction_gain: float = 16.0

    mapping_memory_shear_minimum_ratio: float = 1.0
    mapping_memory_shear_activation: float = 3.0
    mapping_corrector_iterations: int = 5
    mapping_corrector_relaxation: float = 1.0
    mapping_outer_tolerance: float = 1.0e-10

    def __post_init__(self) -> None:
        super().__post_init__()
        if not 0.0 <= self.mapping_bias_onset < self.mapping_bias_full:
            raise ValueError("invalid mapping-backstress bias interval")
        if not 0.0 <= self.mapping_density_full < self.mapping_density_cutoff <= 1.0:
            raise ValueError("invalid mapping-backstress density interval")
        if self.mapping_backstress_rate <= 0.0:
            raise ValueError("mapping backstress rate must be positive")
        if not 0.0 < self.mapping_backstress_capacity_fraction <= 1.0:
            raise ValueError("mapping backstress capacity fraction must lie in (0,1]")
        if not 0.0 < self.mapping_center_limit_ratio < 1.0:
            raise ValueError("mapping center limit ratio must lie in (0,1)")
        if self.mapping_core_radius_ratio <= 0.0:
            raise ValueError("mapping core radius ratio must be positive")
        if not 0.0 <= self.mapping_ray_flow_weight < 1.0:
            raise ValueError("mapping ray-flow weight must lie in [0,1)")
        if min(
            self.mapping_fabric_dilation_rate,
            self.mapping_fabric_recovery_rate,
            self.mapping_fabric_saturation,
        ) <= 0.0:
            raise ValueError("mapping fabric controls must be positive")
        if self.mapping_fabric_flow_weight < 0.0:
            raise ValueError("mapping fabric flow weight must be nonnegative")
        if not 0.0 <= self.mapping_directional_ratchet_weight < 0.5:
            raise ValueError("mapping directional ratchet weight must lie in [0,0.5)")
        if not (
            0.0 <= self.mapping_high_bias_contraction_onset
            < self.mapping_high_bias_contraction_full
        ):
            raise ValueError("invalid high-bias contraction interval")
        if self.mapping_high_bias_contraction_gain < 0.0:
            raise ValueError("high-bias contraction gain must be nonnegative")
        if not 0.0 <= self.mapping_memory_shear_minimum_ratio <= 1.0:
            raise ValueError("mapping minimum shear ratio must lie in [0,1]")
        if self.mapping_memory_shear_activation <= 0.0:
            raise ValueError("mapping memory activation must be positive")
        if not 1 <= self.mapping_corrector_iterations <= 6:
            raise ValueError("mapping corrector requires one to six fixed passes")
        if not 0.0 < self.mapping_corrector_relaxation <= 1.0:
            raise ValueError("mapping corrector relaxation must lie in (0,1]")
        if self.mapping_outer_tolerance <= 0.0:
            raise ValueError("mapping outer tolerance must be positive")


@dataclass
class RIVASandMappingBackstressState(RIVASandLooseBiasedShearFlowState):
    """Host-committed directional state for the translated mapping."""

    mapping_anchor: Tensor = field(
        default_factory=lambda: np.zeros((3, 3), dtype=float)
    )
    mapping_backstress: Tensor = field(
        default_factory=lambda: np.zeros((3, 3), dtype=float)
    )
    mapping_directional_fabric: Tensor = field(
        default_factory=lambda: np.zeros((3, 3), dtype=float)
    )
    mapping_gate_value: float = 0.0
    mapping_capacity: float = 0.0
    mapping_kinematic_denominator: float = 0.0
    mapping_shear_modulus_ratio: float = 1.0
    mapping_phase_contraction_scale: float = 1.0
    mapping_outer_residual: float = 0.0
    mapping_stress_corrections: int = 0
    mapping_corrector_passes: int = 0
    mapping_monotone_caps: int = 0

    def copy(self) -> "RIVASandMappingBackstressState":
        values: dict[str, object] = {}
        for item in fields(self):
            value = getattr(self, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        return RIVASandMappingBackstressState(**values)


class RIVASandMappingBackstressModel(RIVASandLooseBiasedShearFlowModel):
    """Fixed-cost directional mapping and stress-correction prototype."""

    parameters: RIVASandMappingBackstressParameters
    state: RIVASandMappingBackstressState | None

    def __init__(
        self, parameters: RIVASandMappingBackstressParameters | None = None
    ):
        super().__init__(parameters or RIVASandMappingBackstressParameters())

    def initialize(
        self, stress: Tensor, void_ratio: float
    ) -> RIVASandMappingBackstressState:
        base = super().initialize(stress, void_ratio)
        values: dict[str, object] = {}
        for item in fields(RIVASandLooseBiasedShearFlowState):
            value = getattr(base, item.name)
            values[item.name] = value.copy() if isinstance(value, np.ndarray) else value
        self.state = RIVASandMappingBackstressState(**values)
        return self.state.copy()

    def _raw_mapping_gate(self, state: RIVASandState) -> float:
        cfg = self.parameters
        if not cfg.mapping_backstress_enabled or not state.cyclic_phase_active:
            return 0.0
        bias = self._smoothstep(
            (state.static_bias_index - cfg.mapping_bias_onset)
            / (cfg.mapping_bias_full - cfg.mapping_bias_onset)
        )
        density = 1.0 - self._smoothstep(
            (self.initial_relative_density(state) - cfg.mapping_density_full)
            / (cfg.mapping_density_cutoff - cfg.mapping_density_full)
        )
        return float(bias * density)

    def mapping_gate(self, state: RIVASandState) -> float:
        if isinstance(state, RIVASandMappingBackstressState) and state.cyclic_phase_active:
            return float(state.mapping_gate_value)
        return self._raw_mapping_gate(state)

    def begin_cyclic_phase(
        self, *, reference_stress: Tensor | None = None
    ) -> RIVASandMappingBackstressState:
        super().begin_cyclic_phase(reference_stress=reference_stress)
        if not isinstance(self.state, RIVASandMappingBackstressState):
            raise TypeError("mapping-backstress activation lost its state type")
        self.state.mapping_anchor = self.state.alpha.copy()
        self.state.mapping_backstress.fill(0.0)
        self.state.mapping_directional_fabric.fill(0.0)
        self.state.mapping_gate_value = self._raw_mapping_gate(self.state)
        self.state.mapping_capacity = 0.0
        self.state.mapping_kinematic_denominator = 0.0
        self.state.mapping_shear_modulus_ratio = 1.0
        self.state.mapping_phase_contraction_scale = 1.0
        self.state.mapping_outer_residual = 0.0
        self.state.mapping_stress_corrections = 0
        self.state.mapping_corrector_passes = 0
        self.state.mapping_monotone_caps = 0
        return self.state.copy()

    def _mapping_capacity(
        self, state: RIVASandMappingBackstressState, mb: float
    ) -> float:
        cfg = self.parameters
        radius = SQRT23 * mb
        center_limit = cfg.mapping_center_limit_ratio * radius
        clearance = max(center_limit - tensor_norm(state.mapping_anchor), 0.0)
        return float(cfg.mapping_backstress_capacity_fraction * clearance)

    def _mapping_intersection(
        self,
        alpha: Tensor,
        center: Tensor,
        mb: float,
        fallback: Tensor,
    ) -> tuple[float, Tensor, Tensor, int]:
        """Map through a finite directional core to the fixed outer cone."""
        cfg = self.parameters
        radius = SQRT23 * mb
        core = cfg.mapping_core_radius_ratio * radius
        raw_ray = np.asarray(alpha, dtype=float) - np.asarray(center, dtype=float)
        branch = _unit(fallback)
        ray = raw_ray.copy()
        if tensor_norm(branch) > 1.0e-14:
            along = float(np.sum(ray * branch))
            if along < core:
                ray += (core - along) * branch
        if tensor_norm(ray) <= 1.0e-14:
            if tensor_norm(branch) <= 1.0e-14:
                zero = np.zeros((3, 3), dtype=float)
                return cfg.beta0, zero, zero, 1
            ray = core * branch

        qa = float(np.sum(ray * ray))
        qb = 2.0 * float(np.sum(center * ray))
        qc = float(np.sum(center * center)) - radius * radius
        discriminant = qb * qb - 4.0 * qa * qc
        if qa <= 1.0e-16 or discriminant < 0.0:
            direction = _unit(ray)
            return cfg.beta0, direction, direction, 1
        root_term = np.sqrt(max(discriminant, 0.0))
        roots = ((-qb + root_term) / (2.0 * qa), (-qb - root_term) / (2.0 * qa))
        positive = [root for root in roots if root >= 0.0 and np.isfinite(root)]
        if not positive:
            direction = _unit(ray)
            return cfg.beta0, direction, direction, 1
        ray_scale = max(positive)
        mapped = center + ray_scale * ray
        distance = tensor_norm(mapped - alpha)
        beta = max(distance / max(tensor_norm(raw_ray), core), 1.0e-6)
        return float(beta), _unit(mapped), _unit(ray), 0

    def _flow_direction(
        self,
        normal: Tensor,
        ray: Tensor,
        fabric: Tensor,
        bias_direction: Tensor | None = None,
    ) -> Tensor:
        cfg = self.parameters
        fabric_transverse = fabric - float(np.sum(fabric * normal)) * normal
        flow = (
            (1.0 - cfg.mapping_ray_flow_weight) * normal
            + cfg.mapping_ray_flow_weight * ray
            + cfg.mapping_fabric_flow_weight * fabric_transverse
        )
        if bias_direction is not None:
            flow += (
                cfg.mapping_directional_ratchet_weight
                * _unit(bias_direction)
            )
        if tensor_norm(flow) <= 1.0e-14:
            return normal.copy()
        flow = _unit(flow)
        if float(np.sum(flow * normal)) <= 0.05:
            flow = _unit(normal + 0.05 * fabric_transverse)
        return flow

    def _backstress_update(
        self, old: Tensor, flow: Tensor, capacity: float, delta_lambda: float
    ) -> tuple[Tensor, Tensor]:
        """Exact Armstrong--Frederick update for a fixed corrector direction."""
        cfg = self.parameters
        if capacity <= 1.0e-14 or delta_lambda <= 0.0:
            rate = cfg.mapping_backstress_rate * (capacity * flow - old)
            return old.copy(), rate
        target = capacity * flow
        decay = np.exp(-cfg.mapping_backstress_rate * delta_lambda)
        end = target + (old - target) * decay
        midpoint = 0.5 * (old + end)
        rate = cfg.mapping_backstress_rate * (capacity * flow - midpoint)
        return end, rate

    def _fabric_update(
        self,
        old: Tensor,
        normal: Tensor,
        dilatancy: float,
        delta_lambda: float,
    ) -> Tensor:
        cfg = self.parameters
        dilation = max(-dilatancy, 0.0)
        contraction = max(dilatancy, 0.0)
        rate = (
            cfg.mapping_fabric_dilation_rate * dilation
            + cfg.mapping_fabric_recovery_rate * contraction
        )
        if rate <= 1.0e-14 or delta_lambda <= 0.0:
            return old.copy()
        target = (
            -cfg.mapping_fabric_dilation_rate
            * dilation
            * cfg.mapping_fabric_saturation
            * normal
            / rate
        )
        return target + (old - target) * np.exp(-rate * delta_lambda)

    def _memory_shear_ratio(self, backstress: Tensor, capacity: float) -> float:
        cfg = self.parameters
        ratio = tensor_norm(backstress) / max(capacity, 1.0e-14)
        activity = 1.0 - np.exp(-cfg.mapping_memory_shear_activation * ratio * ratio)
        return float(
            1.0 - (1.0 - cfg.mapping_memory_shear_minimum_ratio) * activity
        )

    def _bias_ratchet_increment(
        self, old: RIVASandState, trial: RIVASandState
    ) -> tuple[float, Tensor]:
        increment, direction = super()._bias_ratchet_increment(old, trial)
        gate = self.mapping_gate(old)
        return float((1.0 - gate) * increment), direction

    def _directional_phase_contraction_scale(
        self,
        state: RIVASandMappingBackstressState,
        backstress: Tensor | None = None,
        capacity: float | None = None,
    ) -> float:
        """Reduce high-bias contraction through evolved directional memory."""
        cfg = self.parameters
        high_bias = self._smoothstep(
            (
                state.static_bias_index
                - cfg.mapping_high_bias_contraction_onset
            )
            / (
                cfg.mapping_high_bias_contraction_full
                - cfg.mapping_high_bias_contraction_onset
            )
        )
        direction = _unit(state.static_bias_tensor)
        active_backstress = (
            state.mapping_backstress if backstress is None else backstress
        )
        active_capacity = state.mapping_capacity if capacity is None else capacity
        memory = abs(float(np.sum(active_backstress * direction))) / max(
            active_capacity, 1.0e-14
        )
        return float(np.exp(
            -cfg.mapping_high_bias_contraction_gain
            * self.mapping_gate(state)
            * high_bias
            * memory
        ))

    def _apply_host_phase_volume(self, old, state, deps):
        """Couple PT contraction to the committed directional backstress.

        The inherited continuous pressure wave is retained.  Only the new
        irreversible host increment is scaled, using the average tensorial
        memory at the ends of the host step, and pressure compatibility is
        then rebuilt from the corrected irreversible volume.
        """
        result = super()._apply_host_phase_volume(old, state, deps)
        if (
            not isinstance(old, RIVASandMappingBackstressState)
            or not isinstance(result, RIVASandMappingBackstressState)
            or self.mapping_gate(old) <= 1.0e-14
        ):
            return result
        scale = 0.5 * (
            self._directional_phase_contraction_scale(old)
            + self._directional_phase_contraction_scale(result)
        )
        increment = result.phase_irreversible_volume - old.phase_irreversible_volume
        correction = (scale - 1.0) * increment
        result.phase_irreversible_volume = float(
            old.phase_irreversible_volume + scale * increment
        )
        result.eps_v_irreversible = float(
            result.eps_v_irreversible + correction
        )
        result.mapping_phase_contraction_scale = float(scale)
        self._rebuild_phase_compatible_pressure(result)
        return result

    def _backbone_forward_euler(
        self,
        old: RIVASandState,
        deps: Tensor,
        *,
        force_reversal: bool,
        allow_legacy_reversal: bool,
    ) -> RIVASandState:
        if (
            not isinstance(old, RIVASandMappingBackstressState)
            or self.mapping_gate(old) <= 1.0e-14
            or not old.cyclic_phase_active
        ):
            return super()._backbone_forward_euler(
                old,
                deps,
                force_reversal=force_reversal,
                allow_legacy_reversal=allow_legacy_reversal,
            )

        cfg = self.parameters
        gate = self.mapping_gate(old)
        state = old.copy()
        deps = np.asarray(deps, dtype=float)
        deps_v = float(np.trace(deps))
        deps_d = deps - deps_v * I3 / 3.0
        s_old, p_raw, _ = invariants(old.stress)
        p_old = max(p_raw, cfg.p_min)
        shear_maximum, bulk = self.moduli(p_old)
        mb_old, _, _ = self.surfaces(p_old, old.void_ratio)
        capacity = self._mapping_capacity(old, mb_old)
        shear_ratio = self._memory_shear_ratio(old.mapping_backstress, capacity)
        shear = shear_maximum * shear_ratio
        threshold = self.confining_strain_at_pressure(cfg.p_min, old.pressure_anchor)
        floor_active = cfg.compatibility_enabled and old.eps_v_confining >= threshold
        coupling_bulk = 0.0 if floor_active else bulk

        alpha0 = old.alpha0.copy()
        alpha01 = old.alpha01.copy()
        ep_eq = old.ep_eq_since_reversal
        reversals = old.reversals
        if cfg.compatibility_enabled:
            p_elastic, _ = self.pressure_from_confining_strain(
                old.eps_v_confining + deps_v, old.pressure_anchor
            )
        else:
            p_elastic = max(p_old - bulk * deps_v, cfg.p_min)
        alpha_trial = (s_old + 2.0 * shear * deps_d) / p_elastic
        legacy_measure = float(np.sum((alpha_trial - alpha0) * (alpha_trial - old.alpha)))
        reversal_event = force_reversal or (
            allow_legacy_reversal and legacy_measure < 0.0
        )
        if reversal_event:
            alpha0, alpha01, ep_eq, _, reversals = self._register_reversal(
                old, alpha0, alpha01, ep_eq
            )

        loading_direction = _unit(deps_d)
        center_anchor = (1.0 - gate) * alpha0 + gate * old.mapping_anchor
        backstress_end = old.mapping_backstress.copy()
        fabric_end = old.mapping_directional_fabric.copy()
        delta_lambda = 0.0
        d_ir_mid = max(old.D_ir, 0.0)
        d_re_mid = old.D_re
        flow = old.n.copy()
        if tensor_norm(flow) <= 1.0e-14:
            flow = loading_direction.copy()
        normal = flow.copy()
        beta = max(old.beta, 1.0e-6)
        kinematic_denominator = 0.0
        denominator_floor_hits = old.denominator_floor_hits
        mapping_fallbacks = 0
        monotone_caps = old.mapping_monotone_caps
        at_floor = False

        for _ in range(cfg.mapping_corrector_iterations):
            backstress_mid = 0.5 * (old.mapping_backstress + backstress_end)
            fabric_mid = 0.5 * (old.mapping_directional_fabric + fabric_end)
            total_d = d_ir_mid + d_re_mid
            eps_vc_mid = old.eps_v_confining + 0.5 * (
                deps_v + delta_lambda * total_d
            )
            if cfg.compatibility_enabled:
                p_mid, at_floor_mid = self.pressure_from_confining_strain(
                    eps_vc_mid, old.pressure_anchor
                )
                at_floor = at_floor or at_floor_mid
            else:
                p_mid = max(
                    p_old - 0.5 * bulk * (deps_v + delta_lambda * total_d),
                    cfg.p_min,
                )
                at_floor = at_floor or p_mid <= cfg.p_min
            s_mid = s_old + shear * (deps_d - delta_lambda * flow)
            alpha_mid = s_mid / p_mid
            void_ratio_mid = old.void_ratio + 0.5 * (1.0 + old.void_ratio) * deps_v
            mb_mid, _, _ = self.surfaces(p_mid, void_ratio_mid)
            center_mid = center_anchor + gate * backstress_mid
            beta, normal, ray, fallback = self._mapping_intersection(
                alpha_mid, center_mid, mb_mid, loading_direction
            )
            mapping_fallbacks += fallback
            flow = self._flow_direction(
                normal, ray, fabric_mid, old.static_bias_tensor
            )

            dilatancy_fabric = (
                old.fabric
                + gate * cfg.mapping_fabric_dilatancy_weight * fabric_mid
            )
            d_ir_mid, d_re_mid = self._dilatancy_components(
                alpha_mid,
                flow,
                center_mid,
                beta,
                dilatancy_fabric,
                p_mid,
                void_ratio_mid,
                old.eps_v_irreversible,
                old.eps_v_reversible,
            )
            d_ir_mid *= self.irreversible_contraction_factor(old)
            d_ir_mid *= self._directional_phase_contraction_scale(
                old, backstress_mid, capacity
            )
            if delta_lambda > 0.0 and d_re_mid > 0.0:
                d_re_mid = min(d_re_mid, old.eps_v_reversible / delta_lambda)

            backstress_candidate, backstress_rate = self._backstress_update(
                old.mapping_backstress, flow, capacity, delta_lambda
            )
            normal_flow = max(float(np.sum(normal * flow)), 0.05)
            alpha_normal = float(np.sum(alpha_mid * normal))
            hardening = (
                p_mid
                * self.hardening_prefactor(p_mid)
                * max(beta, 1.0e-12) ** cfg.m
            )
            kinematic_denominator = (
                gate * p_mid * float(np.sum(normal * backstress_rate))
            )
            denominator = (
                2.0 * shear * normal_flow
                + FAC23 * hardening
                + kinematic_denominator
                - coupling_bulk * (d_ir_mid + d_re_mid) * alpha_normal
            )
            denominator_floor = cfg.denominator_floor_ratio * 2.0 * shear
            if denominator < denominator_floor:
                denominator = denominator_floor
                denominator_floor_hits += 1
            numerator = (
                2.0 * shear * float(np.sum(deps_d * normal))
                + coupling_bulk * deps_v * alpha_normal
            )
            candidate = max(0.0, numerator / denominator)
            candidate = (
                (1.0 - cfg.mapping_corrector_relaxation) * delta_lambda
                + cfg.mapping_corrector_relaxation * candidate
            )
            flow_projection = float(np.sum(deps_d * flow))
            if flow_projection > 1.0e-16:
                monotone_limit = float(np.sum(deps_d * deps_d)) / flow_projection
                if candidate > monotone_limit:
                    candidate = monotone_limit
                    monotone_caps += 1
            delta_lambda = float(candidate)
            backstress_end, _ = self._backstress_update(
                old.mapping_backstress, flow, capacity, delta_lambda
            )
            fabric_end = self._fabric_update(
                old.mapping_directional_fabric,
                normal,
                d_ir_mid + d_re_mid,
                delta_lambda,
            )

        s_trial_end = s_old + 2.0 * shear * deps_d
        s_end = s_trial_end - 2.0 * shear * delta_lambda * flow
        eps_v_irreversible = old.eps_v_irreversible - delta_lambda * max(d_ir_mid, 0.0)
        if delta_lambda > 0.0 and d_re_mid > 0.0:
            d_re_mid = min(d_re_mid, old.eps_v_reversible / delta_lambda)
        eps_v_reversible = max(old.eps_v_reversible - delta_lambda * d_re_mid, 0.0)
        eps_v_total = old.eps_v_total + deps_v
        eps_v_confining = eps_v_total - eps_v_irreversible - eps_v_reversible
        if cfg.compatibility_enabled:
            pressure, at_floor_end = self.pressure_from_confining_strain(
                eps_v_confining, old.pressure_anchor
            )
            at_floor = at_floor or at_floor_end
        else:
            pressure = max(
                p_old - bulk * (deps_v + delta_lambda * (d_ir_mid + d_re_mid)),
                cfg.p_min,
            )
            at_floor = at_floor or pressure <= cfg.p_min
        void_ratio = old.void_ratio + (1.0 + old.void_ratio) * deps_v
        mb, _, _ = self.surfaces(pressure, void_ratio)
        radius_stress = pressure * SQRT23 * mb
        alpha = s_end / pressure
        outer_residual = max(tensor_norm(alpha) - SQRT23 * mb, 0.0)
        corrections = 0

        if outer_residual > cfg.mapping_outer_tolerance:
            plastic_vector = 2.0 * shear * flow
            qa = float(np.sum(plastic_vector * plastic_vector))
            qb = -2.0 * float(np.sum(s_trial_end * plastic_vector))
            qc = float(np.sum(s_trial_end * s_trial_end)) - radius_stress * radius_stress
            discriminant = qb * qb - 4.0 * qa * qc
            admissible: list[float] = []
            if qa > 1.0e-16 and discriminant >= 0.0:
                root_term = np.sqrt(max(discriminant, 0.0))
                roots = ((-qb - root_term) / (2.0 * qa), (-qb + root_term) / (2.0 * qa))
                admissible = [root for root in roots if root >= 0.0 and np.isfinite(root)]
            if admissible:
                corrected = min(admissible, key=lambda value: abs(value - delta_lambda))
                flow_projection = float(np.sum(deps_d * flow))
                if flow_projection > 1.0e-16:
                    corrected = min(
                        corrected,
                        float(np.sum(deps_d * deps_d)) / flow_projection,
                    )
                if abs(corrected - delta_lambda) > 1.0e-14:
                    delta_lambda = float(corrected)
                    corrections = 1
                    backstress_end, _ = self._backstress_update(
                        old.mapping_backstress, flow, capacity, delta_lambda
                    )
                    fabric_end = self._fabric_update(
                        old.mapping_directional_fabric,
                        normal,
                        d_ir_mid + d_re_mid,
                        delta_lambda,
                    )
                    s_end = s_trial_end - 2.0 * shear * delta_lambda * flow
                    eps_v_irreversible = (
                        old.eps_v_irreversible
                        - delta_lambda * max(d_ir_mid, 0.0)
                    )
                    if d_re_mid > 0.0:
                        d_re_mid = min(
                            d_re_mid, old.eps_v_reversible / delta_lambda
                        )
                    eps_v_reversible = max(
                        old.eps_v_reversible - delta_lambda * d_re_mid, 0.0
                    )
                    eps_v_confining = (
                        eps_v_total - eps_v_irreversible - eps_v_reversible
                    )
                    if cfg.compatibility_enabled:
                        pressure, at_floor_return = self.pressure_from_confining_strain(
                            eps_v_confining, old.pressure_anchor
                        )
                        at_floor = at_floor or at_floor_return
                    else:
                        pressure = max(
                            p_old
                            - bulk * (
                                deps_v + delta_lambda * (d_ir_mid + d_re_mid)
                            ),
                            cfg.p_min,
                        )
                    alpha = s_end / pressure
                    mb, _, _ = self.surfaces(pressure, void_ratio)
                    outer_residual = max(tensor_norm(alpha) - SQRT23 * mb, 0.0)

        stress = s_end - pressure * I3
        center_end = center_anchor + gate * backstress_end
        beta_end, normal_end, ray_end, fallback = self._mapping_intersection(
            alpha, center_end, mb, flow
        )
        mapping_fallbacks += fallback
        flow_end = self._flow_direction(
            normal_end, ray_end, fabric_end, old.static_bias_tensor
        )
        dilatancy_fabric_end = (
            old.fabric + gate * cfg.mapping_fabric_dilatancy_weight * fabric_end
        )
        d_ir_end, d_re_end = self._dilatancy_components(
            alpha,
            flow_end,
            center_end,
            beta_end,
            dilatancy_fabric_end,
            pressure,
            void_ratio,
            eps_v_irreversible,
            eps_v_reversible,
        )
        d_ir_end *= self.irreversible_contraction_factor(old)
        d_ir_end *= self._directional_phase_contraction_scale(
            old, backstress_end, capacity
        )

        ep_eq += delta_lambda
        state.stress = stress
        state.alpha = alpha
        state.alpha0 = alpha0
        state.alpha01 = alpha01
        state.n = flow_end
        state.fabric = old.fabric.copy()
        state.D_ir = float(d_ir_end)
        state.D_re = float(d_re_end)
        state.D = float(d_ir_end + d_re_end)
        state.beta = float(beta_end)
        state.lambda_total = float(old.lambda_total + delta_lambda)
        state.ep_eq_since_reversal = float(ep_eq)
        state.void_ratio = float(void_ratio)
        state.reversals = reversals
        state.pressure_floor_hits = old.pressure_floor_hits + int(at_floor)
        state.denominator_floor_hits = denominator_floor_hits
        state.beta_fallbacks = old.beta_fallbacks + mapping_fallbacks
        state.eps_v_total = float(eps_v_total)
        state.eps_v_confining = float(eps_v_confining)
        state.eps_v_irreversible = float(eps_v_irreversible)
        state.eps_v_reversible = float(eps_v_reversible)
        state.mapping_backstress = backstress_end
        state.mapping_directional_fabric = fabric_end
        state.mapping_capacity = float(capacity)
        state.mapping_kinematic_denominator = float(kinematic_denominator)
        state.mapping_shear_modulus_ratio = float(shear_ratio)
        state.mapping_outer_residual = float(outer_residual)
        state.mapping_stress_corrections = (
            old.mapping_stress_corrections + corrections
        )
        state.mapping_corrector_passes = cfg.mapping_corrector_iterations
        state.mapping_monotone_caps = monotone_caps

        if state.reversals > old.reversals:
            reversal_deviator, _, _ = invariants(old.stress)
            excursion = tensor_norm(reversal_deviator - old.last_reversal_deviator)
            divisor = 1.0 if old.amplitude_reversals == 0 else 2.0
            amplitude = excursion / (divisor * max(old.pressure_anchor, cfg.p_min))
            if amplitude > 1.0e-12:
                state.cyclic_amplitude = float(amplitude)
                state.amplitude_factor = self.amplitude_factor_for_state(amplitude, old)
            state.last_reversal_deviator = reversal_deviator.copy()
            state.amplitude_reversals = old.amplitude_reversals + 1
        return state

    def _host_outer_stress_correction(
        self, state: RIVASandMappingBackstressState
    ) -> RIVASandMappingBackstressState:
        """Return the post-compatibility stress to the outer cone.

        The signed phase-volume mechanism is committed once per host step by
        the parent class.  Its pressure correction can move an otherwise
        admissible deviator outside the pressure-dependent cone.  Complete
        the same explicit update with an isochoric plastic multiplier return;
        this preserves the compatible pressure and avoids radial stress
        clipping without a corresponding plastic-state increment.
        """
        if self.mapping_gate(state) <= 1.0e-14:
            return state
        cfg = self.parameters
        corrected = state.copy()
        correction_count = 0
        for _ in range(2):
            deviator, pressure_raw, _ = invariants(corrected.stress)
            pressure = max(pressure_raw, cfg.p_min)
            mb, _, _ = self.surfaces(pressure, corrected.void_ratio)
            radius = pressure * SQRT23 * mb
            norm = tensor_norm(deviator)
            residual = max(norm / pressure - SQRT23 * mb, 0.0)
            corrected.mapping_outer_residual = float(residual)
            if residual <= cfg.mapping_outer_tolerance:
                break
            direction = _unit(deviator)
            shear, _ = self.moduli(pressure)
            delta_lambda = max((norm - radius) / max(2.0 * shear, 1.0e-14), 0.0)
            if delta_lambda <= 0.0:
                break
            corrected.mapping_backstress, _ = self._backstress_update(
                corrected.mapping_backstress,
                direction,
                corrected.mapping_capacity,
                delta_lambda,
            )
            corrected.mapping_directional_fabric = self._fabric_update(
                corrected.mapping_directional_fabric,
                direction,
                corrected.D,
                delta_lambda,
            )
            deviator = radius * direction
            corrected.stress = deviator - pressure * I3
            corrected.alpha = deviator / pressure
            corrected.lambda_total = float(corrected.lambda_total + delta_lambda)
            corrected.ep_eq_since_reversal = float(
                corrected.ep_eq_since_reversal + delta_lambda
            )
            center = (
                (1.0 - self.mapping_gate(corrected)) * corrected.alpha0
                + self.mapping_gate(corrected)
                * (corrected.mapping_anchor + corrected.mapping_backstress)
            )
            beta, normal, ray, fallback = self._mapping_intersection(
                corrected.alpha, center, mb, direction
            )
            corrected.beta = float(beta)
            corrected.n = self._flow_direction(
                normal,
                ray,
                corrected.mapping_directional_fabric,
                corrected.static_bias_tensor,
            )
            corrected.beta_fallbacks += fallback
            correction_count += 1
        corrected.mapping_stress_corrections += correction_count
        deviator, pressure_raw, _ = invariants(corrected.stress)
        pressure = max(pressure_raw, cfg.p_min)
        mb, _, _ = self.surfaces(pressure, corrected.void_ratio)
        corrected.mapping_outer_residual = float(max(
            tensor_norm(deviator) / pressure - SQRT23 * mb, 0.0
        ))
        return corrected

    def advance_fixed(
        self,
        initial: RIVASandMappingBackstressState,
        deps: Tensor,
        nsub: int = 1,
    ):
        state, info = super().advance_fixed(initial, deps, nsub)
        if not isinstance(state, RIVASandMappingBackstressState):
            raise TypeError("mapping-backstress update lost its state type")
        if self.mapping_gate(state) > 1.0e-14:
            state = self._host_outer_stress_correction(state)
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
        if isinstance(current, RIVASandMappingBackstressState):
            values.update(
                mapping_gate=self.mapping_gate(current),
                mapping_backstress_norm=tensor_norm(current.mapping_backstress),
                mapping_fabric_norm=tensor_norm(current.mapping_directional_fabric),
                mapping_capacity=current.mapping_capacity,
                mapping_kinematic_denominator=current.mapping_kinematic_denominator,
                mapping_shear_modulus_ratio=current.mapping_shear_modulus_ratio,
                mapping_phase_contraction_scale=(
                    current.mapping_phase_contraction_scale
                ),
                mapping_outer_residual=current.mapping_outer_residual,
                mapping_stress_corrections=float(current.mapping_stress_corrections),
                mapping_corrector_passes=float(current.mapping_corrector_passes),
                mapping_monotone_caps=float(current.mapping_monotone_caps),
            )
        return values


__all__ = [
    "RIVASandMappingBackstressModel",
    "RIVASandMappingBackstressParameters",
    "RIVASandMappingBackstressState",
]
