# bias_reversible_volume under static bias at low amplitude / low confinement

Found 2026-08-24/25 while running the LEAP-Asia-2019 sloping-ground slice BVP
(RPI-A, 125x25x1 SSPbrickUP, prototype scale). Reported by Yu-Wei Hwang /
Claude. This note documents the defect, the contained fix in this PR, its
validation, and a proposed V9 reformulation.

## Symptom

A 5-degree slope model with an admissible initial state (K0 = 0.485 +
infinite-slope shear, severity <= 0.84 everywhere) collapses within the first
seconds of a ramped 1 Hz motion while the input is still at micro amplitude
(0.0002 g): sigma'_v breathes with a growing envelope across the biased
region, the solver cascades, and the slope "statically liquefies". A
level-ground twin (no static bias) runs the full motion flawlessly.

Minimal single-element reproduction (`tests/RIVASand_bias_micro_cycle.tcl`,
run it against a pre-patch build to see the failure):

- sigma_v0 = 25 kPa, alpha = tau_s/sigma_v0 = 0.09, tau_cyc = 0.005 kPa
  (CSR 2e-4), constant-volume DSS: sigma'_v oscillates 25 -> 21.8/31.5 ->
  ... -> 13.5/47.8 kPa within five cycles while the shear strain stays frozen
  at ~3.5e-5 and eps_v is constant. alpha = 0: perfectly quiet. Independent
  of steps/cycle (20/200/2000) and nSub (1/8).
- State trace: `bias_reversible_volume` (state[86]) and `eps_v_confining`
  (state[46]) ratchet each reversal; envelope follows 1-exp(-reversals/6).

## Root cause (`RIVASandKernel.h`, `riva_bias_reversible_volume_target`)

1. `phase = clip(dynamic_projection/(pressure_anchor*cyclic_amplitude))` is
   NORMALIZED by the cyclic amplitude, so the oscillation term swings its
   full calibrated magnitude at ANY amplitude — the mechanism has no
   amplitude consistency. (Contrast: PM4Sand and SANISAND drive all cyclic
   volumetric machinery by the plastic multiplier — dFabric ~ <L><-D> in
   PM4Sand.cpp:1338-1342, h = b0/|(alpha-alpha_in):n| in
   ManzariDafalias.cpp:2929 — so their response vanishes with vanishing
   amplitude by construction.)
2. `bias_factor = (0.646/bias)^2` DIVERGES as the projected bias decreases;
   gentle-slope states (bias ~ 0.19, below the calibrated envelope's floor of
   0.229) are amplified ~11x.
3. At low pressure anchors the fixed volumetric amplitude (7.6e-5-scale)
   converts through the small bulk modulus into confinement-scale pressure
   swings (sigma_v0 = 2.5 kPa, CSR 0.05: +/-50% sigma'_v per cycle).

## Fix in this PR (calibration-preserving containment)

Two gates in `riva_bias_reversible_volume_target`, both exactly 1.0 on every
calibrated loading path so the frozen V8 oracle is preserved bit-for-bit:

- amplitude gate: `smoothstep((amp*max(anchor,pmin)/max(p,pmin) - 0.12)/
  (0.35 - 0.12))`. The discriminant is normalized by the CURRENT pressure so
  post-liquefaction golden cycles (legitimately collapsed stress amplitude,
  collapsed p') stay fully open; the golden floor of the discriminant is
  0.5375, the window was chosen empirically as the widest that keeps all four
  bias/mixed golden cases bit-identical (0.15/0.40 already breaks
  mixed_3d_bias015 through provisional-state excursions).
- low-pressure fade: `smoothstep(anchor/12 kPa)` (golden anchors are all
  26.27 kPa).

Plus an adapter escape hatch `-noBiasVolume` (RIVASand.cpp/.h) that disables
the block entirely; element-level A/B at calibrated amplitudes shows N-to-3%
and strain accumulation identical with it off — only peak r_u oscillation
amplitudes differ.

Validation: GoldenReplay 4/4 bias+mixed cases bit-PASS with the gates; the
`cyclic_zero_bias_reference` case shows a pre-existing 1.8087e-6 last-step
tolerance miss that is IDENTICAL with and without this change (present on the
unpatched 2026-08-24 build as well — likely a compiler/flag drift since the
oracle freeze, worth a separate look). KernelStateTest passes. The
new regression `tests/RIVASand_bias_micro_cycle.tcl` passes (pre-patch it
fails with ~90% sigma'_v drift).

## Remaining hole and V9 proposal

The gates do not cover the mid-band (discriminant 0.12-0.5) where slope
states at bias ~0.19 still receive a near-full-size oscillation, and the
p-normalized discriminant forms a positive feedback (p drops -> gate opens ->
mechanism drops p further) — observed as premature full-depth ru in the BVP
at 0.02 g. Because the golden envelope genuinely overlaps this band (bias
floor 0.229 vs slope 0.19), no gate variable separates them robustly.

Proposed V9 reformulation (PM4Sand/SANISAND-consistent): scale the target by
the plastic activity of the cycle, e.g. multiply by
`smoothstep(ep_half/ep_ref)` where `ep_half` is the plastic multiplier
accumulated over the last completed half-cycle (new state variable captured
from `ep_eq_since_reversal` at reversal registration; ep_ref calibrated from
the golden histories so calibrated paths sit at 1). Amplitude consistency
then holds structurally and both gates in this PR can be removed. Costs:
state vector 93 -> 94, state_schema/adapter/state-test updates, golden
re-freeze (V9 oracle), Hercules port sync.

## V9 experiment (2026-08-25): plastic-activity gate — result and a decisive negative

Implemented `ep_half_last` (new state field, plastic multiplier of the last
completed half-cycle, captured at reversal registration; state vector 93->94)
and gated the target by `smoothstep(ep_half_last/ep_ref)`, ep_ref = 1.5e-5,
below the golden floor (3.4e-5, weakest half-cycle of cyclic_bias0375_dense).
Golden replay stays bit-identical (the gate saturates at 1.0 on every
calibrated evaluation) — no oracle re-freeze needed.

Findings:
- The plastic gate alone does NOT stabilize sloping-ground BVPs: gravity
  redistribution on a slope is genuinely plastic (ep_half >> 1.5e-5) at tiny
  STRESS amplitude, so the mechanism re-ignites during the gravity settle.
- The stress-amplitude gate alone leaves the p-collapse feedback (p drops ->
  amp*anchor/p rises -> mechanism opens -> p drops further).
- Both gates TOGETHER (plus the low-pressure fade) are strictly safer and all
  evaluate to exactly 1.0 on the golden paths; this is what the branch now
  carries. BVP work still runs `-noBiasVolume`.

The structural conclusion: a SATURATING gate can never be simultaneously
bit-preserving (=1 across the golden envelope) and proportional below it,
because the golden envelope itself contains a weak-plasticity half-cycle
(ep = 3.4e-5) that the frozen calibration treats at full mechanism strength.
Making the response truly amplitude-consistent (PM4Sand-style, response
proportional to the cycle plastic multiplier, e.g. min(1, ep_half/ep_cal)
with ep_cal ~ 3e-3) necessarily rescales that golden half-cycle — i.e. the
proper fix REQUIRES re-freezing the oracle and re-checking the DSS
calibration. That decision belongs to the model owner; the infrastructure
for it (ep_half_last state + gate plumbing) is in place on this branch.
