"""Minimal biased-loop research extension of the frozen RIVA-Sand oracle.

The production equations are unchanged unless ``branch_compliance_enabled``
is active.  The extension applies a reversal-anchored Masing-type tangent
multiplier only to dense, statically biased cyclic states.  It is deliberately
confined to shear compliance: plastic hardening, dilatancy, pressure
compatibility, pore-pressure generation, and the calibrated ratchet are the
frozen RIVA-Sand implementation.

This is private research code, not a production RIVA-Sand release.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rivasand_port.model import (
    RIVASandModel,
    RIVASandParameters,
    RIVASandState,
    invariants,
    tensor_norm,
)


@dataclass(frozen=True)
class RIVASandBiasedLoopParameters(RIVASandParameters):
    """Frozen RIVA-Sand controls plus a compact biased-loop compliance law."""

    # Small global retune of the existing cyclic compliance.  The production
    # value is 0.85; 0.87 corrects the reference-density loop range and damping
    # without changing the observed eight-case cycle event classifications.
    cyclic_shear_modulus_reduction: float = 0.87
    branch_compliance_enabled: bool = True
    branch_compliance_gain: float = 4.5
    branch_compliance_exponent: float = 1.0
    branch_compliance_bias_reference: float = 0.5384061785684372
    branch_compliance_bias_exponent: float = 1.0
    branch_compliance_minimum: float = 0.10
    branch_directional_balance: float = 0.03
    branch_balance_bias_exponent: float = 0.20

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.branch_compliance_gain < 0.0:
            raise ValueError("branch compliance gain must be nonnegative")
        if self.branch_compliance_exponent <= 0.0:
            raise ValueError("branch compliance exponent must be positive")
        if self.branch_compliance_bias_reference <= 0.0:
            raise ValueError("branch compliance bias reference must be positive")
        if self.branch_compliance_bias_exponent < 0.0:
            raise ValueError("branch compliance bias exponent must be nonnegative")
        if not 0.0 < self.branch_compliance_minimum <= 1.0:
            raise ValueError("branch compliance minimum must lie in (0,1]")
        if not np.isfinite(self.branch_directional_balance):
            raise ValueError("branch directional balance must be finite")
        if self.branch_balance_bias_exponent < 0.0:
            raise ValueError("branch balance bias exponent must be nonnegative")


class RIVASandBiasedLoopModel(RIVASandModel):
    """RIVA-Sand with one smooth, reversal-anchored shear-compliance modifier."""

    parameters: RIVASandBiasedLoopParameters

    def __init__(self, parameters: RIVASandBiasedLoopParameters | None = None):
        super().__init__(parameters or RIVASandBiasedLoopParameters())

    def branch_progress(self, state: RIVASandState) -> float:
        """Normalized deviatoric excursion from the current reversal point."""
        cfg = self.parameters
        if (
            not cfg.branch_compliance_enabled
            or not state.cyclic_phase_active
            or state.amplitude_reversals < 1
            or state.cyclic_amplitude <= 1.0e-14
        ):
            return 0.0
        deviator, _, _ = invariants(state.stress)
        excursion = tensor_norm(deviator - state.last_reversal_deviator)
        half_cycle_span = (
            2.0
            * max(state.pressure_anchor, cfg.p_min)
            * state.cyclic_amplitude
        )
        return float(np.clip(excursion / max(half_cycle_span, 1.0e-14), 0.0, 1.0))

    def branch_compliance_multiplier(self, state: RIVASandState) -> float:
        """Return 1 at reversal and soften smoothly across each half-cycle."""
        cfg = self.parameters
        if not cfg.branch_compliance_enabled:
            return 1.0
        density = self.dense_state_weight(state)
        bias = self.projected_bias(state)
        if density <= 1.0e-14 or bias <= 1.0e-14:
            return 1.0
        bias_weight = (
            bias / cfg.branch_compliance_bias_reference
        ) ** cfg.branch_compliance_bias_exponent
        progress = self.branch_progress(state)
        direction_norm = tensor_norm(state.cyclic_direction)
        signed_projection = float(
            np.sum(state.static_bias_tensor * state.cyclic_direction)
            / max(direction_norm, 1.0e-14)
        )
        # A weak power of bias lets the stronger alpha=0.375 histories receive
        # more balancing than alpha=0.25 without the severe over-correction
        # produced by a linear bias multiplier.
        signed_bias = float(
            np.sign(signed_projection)
            * min(
                (
                    abs(signed_projection)
                    / cfg.branch_compliance_bias_reference
                ) ** cfg.branch_balance_bias_exponent,
                1.5,
            )
        )
        # Loading in the static-bias direction is slightly stiffer than the
        # reverse branch.  The exponential is smooth, strictly positive, and
        # equals one when the optional balancing control is zero.
        # ``cyclic_direction`` stores the completed excursion used to anchor
        # the current branch, so its sign is opposite to the active loading
        # direction.  A positive stored projection therefore identifies the
        # reverse branch and receives the larger compliance.
        branch_balance = np.exp(cfg.branch_directional_balance * signed_bias)
        compliance = (
            cfg.branch_compliance_gain
            * density
            * bias_weight
            * progress**cfg.branch_compliance_exponent
            * branch_balance
        )
        return float(np.clip(
            1.0 / (1.0 + compliance),
            cfg.branch_compliance_minimum,
            1.0,
        ))

    def moduli_for_state(
        self, pressure: float, state: RIVASandState
    ) -> tuple[float, float]:
        shear, bulk = super().moduli_for_state(pressure, state)
        shear *= self.branch_compliance_multiplier(state)
        return float(shear), float(bulk)

    def dss_history_values(
        self, state: RIVASandState | None = None
    ) -> dict[str, float]:
        current = state or self.state
        if current is None:
            raise RuntimeError("initialize the model first")
        values = super().dss_history_values(current)
        values.update(
            branch_progress=self.branch_progress(current),
            branch_compliance_multiplier=self.branch_compliance_multiplier(current),
            branch_bias_sign=float(np.sign(
                np.sum(current.static_bias_tensor * current.cyclic_direction)
            )),
        )
        return values


__all__ = [
    "RIVASandBiasedLoopModel",
    "RIVASandBiasedLoopParameters",
]
