# Restrained Phase-Transformation Research Checkpoint

Status: private research successor; not production RIVA-Sand and not approved
for OpenSees or Hercules.

Ancestry: exact branch base `origin/research/rivasand-geostatic` at
`51d93d884`. The implementation starts directly from the frozen RIVA-Sand
oracle. Setting `phase_transformation_enabled=False` reproduces every
production state field exactly.

## Purpose

Production RIVA-Sand gives narrow, nearly elastic ribbons in the dense biased
PRJ-4666 DSS tests and a folded, hourglass-like effective-stress path in the
intermediate-density `3484_b025` test. This checkpoint introduces one
restrained phase-transformation mechanism to couple:

1. low-signed-shear compliance and S-shaped loop curvature;
2. dense biased plastic hardening;
3. compatible irreversible contraction;
4. a continuously evolved reversible volume and pressure wave;
5. smooth accumulated phase memory at high static bias; and
6. a density-dependent fixed-direction pressure wave extending from
   intermediate to dense sand.

The model does not assign pore pressure or `ru`. It modifies compatible
volumetric strain, from which effective pressure is rebuilt using the frozen
RIVA-Sand pressure--confining-strain relation.

## Activation and phase coordinate

Let `A` be the normalized, fixed tensor direction of the admitted static-bias
stress and `alpha=s/p` the current stress-ratio tensor. The signed coordinate
and transformation-zone weight are

```
eta_A = alpha : A
eta_PT = phase_ratio * sqrt(2/3) * M_d
chi = tanh(eta_A / (phase_width * eta_PT))
z_PT = max(1 - chi^2, 0)
```

Two activation gates are deliberately separated. The dense shear/hardening
gate remains the product of the frozen dense-state weight and a cubic
smoothstep in projected static bias. Thus, the stress--strain calibration is
unchanged for intermediate-density sand.

The compatible PT-volume gate uses a new density weight

```
w_D = smoothstep((D_r - 0.53) / (0.90 - 0.53))
g_v = w_D * smoothstep(bias / bias_onset)
```

It is nonzero at the `3484_b025` density (`D_r=0.663`, `w_D=0.295`) and reaches
the existing dense calibration at `D_r=0.90`. A separate replacement gate
transitions from zero at `D_r=0.53` to one at `D_r=0.66`; it continuously
removes the frozen branch-relative reversible wave before that wave folds into
an hourglass. Unbiased states and loose states below the onset still execute
the production update exactly.

## Shear response

For each accepted reversal, branch progress is the deviatoric stress excursion
from the last reversal normalized by twice the detected half-cycle amplitude.
A smooth monotonic transition in branch progress is combined with the fixed
phase-zone bell:

```
c_PT = density_weight * bias_weight *
       (phase_compliance_peak * branch_transition * branch_balance
        + phase_compliance_bell_gain * z_PT)
G_PT = G_RIVA / (1 + c_PT)
```

The multiplier is bounded below by `branch_compliance_minimum`. The bell gives
the low-stress portion of each branch its S-shape; the monotonic term preserves
finite loop area after reversal.

The dense biased cyclic shear cut is blended from the frozen value only inside
the same activation gate. Outside the gate the production modulus is exact.

## Restrained plastic hardening and phase memory

Inside the dense biased branch, the static-bias intercept is expressed as

```
H_bias,PT = 87.291302 * bias^3.2 * crossing_decay
            + frozen_margin_term
```

This leaves the alpha=0.25 calibration essentially unchanged while increasing
relative restraint at alpha=0.375. At high bias only, accumulated compatible
phase contraction evolves the smooth memory

```
W_PT = 1 - exp(-abs(eps_v,PT^ir) / phase_memory_reference_volume)
H_bias = H_bias,PT * exp(-phase_memory_hardening_reduction
                        * density_weight * z_PT
                        * high_bias_gate * W_PT)
```

This is not a cycle counter and has no instantaneous amplitude knee. It fades
excess high-bias hardening continuously as actual irreversible phase volume is
accumulated, allowing delayed cyclic mobility while retaining the cycle-11
loop shape.

## Compatible volumetric update

For one host strain increment, phase memory is advanced once after the fixed
constitutive substeps. Let `dell` be the deviatoric strain-path length,
`a_ir` the irreversible density/bias/amplitude/pressure activity, and `s_PT`
the smooth contraction--dilation sign from the mapped loading stress ratio.
Then

```
d eps_v,PT^ir = -phase_contraction_rate * a_ir * z_PT
                * 0.5 * (1 - s_PT) * dell
```

The reversible wave and its accumulated mean shift use separate activities.
This prevents pressure-wave amplitude, mean pressure, and phase lag from being
forced through one scalar multiplier. Let `b_D` be a smooth transition that is
zero at the `3484_b025` reference density (`D_r=0.663`) and one at `D_r=0.90`.
The intermediate-density controls are blended continuously into the dense
calibration:

```
k_wave = 1.1333 + b_D * (1 - 1.1333)
k_mean = 10.7692 + b_D * (1 - 10.7692)
k_relax = 0.10417 + b_D * (1 - 0.10417)

T_i = phase_reversible_scale * k_wave * a_wave,i
      * (chi_i - chi_anchor)
      + phase_reversible_mean_scale * k_mean * a_mean,i * W_PT,i

T_bar = 0.5 * (T_n + T_n+1)

L_relax = phase_reversible_relaxation_strain * k_relax
          * bias_ratio^phase_reversible_relaxation_bias_exponent
r = 1 - exp(-dell / L_relax)
eps_v,PT,n+1^re = eps_v,PT,n^re
                  + r * (T_bar - eps_v,PT,n^re)
```

The relaxation length is multiplied by a smooth static-bias power. The wave,
mean, and irreversible activities also have independently calibrated pressure
and bias powers. The trapezoidal endpoint target is important: evaluating
only the new endpoint produced an artificial timestep dependence in the
effective-stress-path area.

The new phase volumes enter the existing reversible--irreversible
compatibility split. Effective pressure is then reconstructed, and the
deviatoric stress is retained. Internal constitutive substeps see the committed
host-level phase volume, so phase memory is not multiplied by `nSub`. The
intermediate multipliers affect compatible reversible volume only; the dense
shear compliance, hardening, and irreversible contraction laws are not turned
on at `D_r=0.663`.

## Calibrated internal controls

| Control | Value |
|---|---:|
| `phase_cyclic_shear_modulus_reduction` | 0.813 |
| `phase_compliance_peak` | 6.0 |
| `phase_compliance_bell_gain` | 3.5 |
| `phase_compliance_shape` | 1.55 |
| `phase_compliance_location` | 0.442 |
| `phase_compliance_half_width` | 0.787 |
| `branch_directional_balance` | 0.050 |
| `branch_balance_bias_exponent` | 4.0 |
| `branch_balance_bias_cap` | 3.9 |
| `phase_bias_hardening_intercept` | 87.291302 |
| `phase_bias_hardening_exponent` | 3.2 |
| `phase_memory_hardening_reduction` | 2.8 |
| `phase_memory_reference_volume` | 5.0e-5 |
| `phase_ratio` | 0.62 |
| `phase_width` | 0.50 |
| `phase_contraction_rate` | 4.0e-4 |
| `phase_reversible_scale` | -3.0e-3 |
| `phase_reversible_mean_scale` | 6.5e-4 |
| `phase_reversible_relaxation_strain` | 0.0048 |
| `phase_reversible_relaxation_bias_exponent` | 1.6 |
| `phase_potential_anchor_fraction` | 1.0 |
| `phase_pressure_exponent` | 1.6 |
| `phase_wave_pressure_exponent` | 1.0 |
| `phase_mean_pressure_exponent` | 0.0 |
| `phase_amplitude_onset_ratio` | 0.55 |
| `phase_amplitude_full_ratio` | 0.90 |
| `phase_volume_density_onset` | 0.53 |
| `phase_volume_density_full` | 0.90 |
| `phase_volume_replacement_density_full` | 0.66 |
| `phase_intermediate_wave_multiplier` | 1.133333 |
| `phase_intermediate_mean_multiplier` | 10.769231 |
| `phase_intermediate_relaxation_ratio` | 0.104167 |
| `phase_bias_exponent` | 2.0 |
| `phase_wave_bias_exponent` | -0.25 |
| `phase_mean_bias_exponent` | 12.0 |

All production RIVA-Sand parameters retain their frozen values. The table is
an internal research calibration, not a proposed user input list.

## Update order

For each host increment:

1. Copy the committed state and detect the host-level reversal once.
2. If both the dense shear/hardening gate and compatible PT-volume gate are
   zero, execute the frozen RIVA-Sand update.
3. Otherwise, hold committed phase volume fixed during constitutive substeps.
4. At each substep, evaluate pressure-dependent moduli, the PT compliance,
   mapped-surface plastic flow, restrained bias hardening, and frozen state
   evolution.
5. After all substeps, evaluate the total host deviatoric path length.
6. Advance irreversible phase volume once.
7. Relax reversible phase volume once toward its signed PT target.
8. Insert both volumes into the frozen compatibility split.
9. Rebuild effective pressure, retain deviatoric stress, and commit all state.

## Validation result

Cycle-11 values at 32 points/cycle and four constitutive substeps:

| Case | Quantity | Experiment | Production | PT checkpoint |
|---|---|---:|---:|---:|
| PRJ-3484 alpha=0.25 | cycles to criterion | 66.23 | 27.16 | 26.22 |
|  | strain center (%) | 3.337 | 3.375 | 3.388 |
|  | strain-range ratio | 1.000 | 0.657 | 0.657 |
|  | loop-area ratio | 1.000 | 0.911 | 0.912 |
|  | damping ratio | 0.212 | 0.294 | 0.294 |
|  | affine shape RMSE | 0.000 | 0.101 | 0.102 |
|  | mean `ru` | 0.406 | 0.445 | 0.394 |
| PRJ-4666 alpha=0.25 | cycles to criterion | 36.23 | 28.25 | 33.25 |
|  | strain center (%) | 2.087 | 3.574 | 2.292 |
|  | strain-range ratio | 1.000 | 0.210 | 0.811 |
|  | loop-area ratio | 1.000 | 0.131 | 0.735 |
|  | damping ratio | 0.139 | 0.091 | 0.132 |
|  | affine shape RMSE | 0.000 | 0.246 | 0.146 |
|  | mean `ru` | 0.464 | 0.212 | 0.436 |
| PRJ-4666 alpha=0.375 | cycles to criterion | 29.19 | 24.22 | 29.22 |
|  | strain center (%) | 2.568 | 4.809 | 2.636 |
|  | strain-range ratio | 1.000 | 0.197 | 0.979 |
|  | loop-area ratio | 1.000 | 0.126 | 0.767 |
|  | damping ratio | 0.127 | 0.081 | 0.099 |
|  | affine shape RMSE | 0.000 | 0.263 | 0.178 |
|  | mean `ru` | 0.214 | -0.079 | 0.204 |

The effective-stress path is now an explicit calibration target:

| Case | Path quantity | Experiment | Production | PT checkpoint |
|---|---|---:|---:|---:|
| PRJ-3484 alpha=0.25 | mean sigma_v_eff (kPa) | 23.75 | 22.19 | 24.24 |
|  | sigma_v_eff range (kPa) | 11.63 | 8.75 | 11.65 |
|  | enclosed path area (kPa^2) | 44.34 | 6.39 | 42.55 |
|  | correlation of tau and sigma_v_eff | 0.947 | -0.111 | 0.927 |
| PRJ-4666 alpha=0.25 | mean sigma_v_eff (kPa) | 21.45 | 31.54 | 22.56 |
|  | sigma_v_eff range (kPa) | 37.26 | 11.01 | 37.35 |
|  | enclosed path area (kPa^2) | 237.87 | 13.26 | 242.89 |
|  | correlation of tau and sigma_v_eff | 0.939 | 0.006 | 0.952 |
| PRJ-4666 alpha=0.375 | mean sigma_v_eff (kPa) | 31.43 | 43.15 | 31.82 |
|  | sigma_v_eff range (kPa) | 48.16 | 6.72 | 49.35 |
|  | enclosed path area (kPa^2) | 420.77 | 25.47 | 425.85 |
|  | correlation of tau and sigma_v_eff | 0.954 | approximately 0 | 0.944 |

The four unbiased/loose histories (`3484_u015`, `3484_u021`,
`4666_loose_b015`, and `4666_dense_u038`) are numerically identical to
production in every common recorded field. The `3484_b030` path activates the
new volume law but still stops at `N=1.0625`, before its cycle-2 comparison.
The aggregate strict-audit score is 2.239 because the intentionally untouched
production limitations remain and two early-failure histories do not reach
their comparison cycles. All eight cases still fail at least one hard gate.
This checkpoint therefore is not a replacement for production.

## Objectivity and restart

The state/restart test passes with exact restart, zero compatibility residual,
finite response, and exact production equivalence when disabled. The generic
2-versus-4-substep maximum stress difference is 8.070%.

For `3484_b025`, changing from 32 points/cycle with four substeps to 64
points/cycle with two substeps changes cycles to criterion by 3.75%, path range
by 1.18%, path area by 0.58%, mean `ru` by 0.65%, peak-to-peak `ru` by 1.18%,
and stress--strain loop area by 0.23%.

For the two dense biased audits, the corresponding changes are 3.01% and
3.42% in cycles to criterion; 2.30% and 0.10% in path range; 0.23% and 5.81%
in path area; 6.04% and 3.71% in mean `ru`; 2.30% and 0.10% in peak-to-peak
`ru`; and 4.41% and 6.74% in stress--strain loop area.

## Cost and limitations

A Python eight-case audit is about 2.17 times the frozen
production model. This includes Python dispatch and is not a native C++/GPU
kernel benchmark, but it is too expensive to call production-ready without a
native implementation and profiling.

Remaining limitations:

- Stress--strain loop area and damping remain low, especially at alpha=0.375,
  although both are greatly improved over production.
- `ru` waveform history RMSE remains outside the strict gate, particularly at
  alpha=0.375.
- The `3484_b025` cycle-11 effective-stress path no longer has the production
  hourglass and closely matches its experimental mean, range, area, and
  orientation. However, liquefaction remains much too fast (`N=26.22` versus
  `66.23`). A scalar contraction reduction can delay liquefaction but destroys
  the already-correct cycle-11 strain center and loop; this requires a separate
  late-cycle shakedown/accumulation law.
- `3484_b030` still stops at `N=1.0625`, so its cycle-2 comparison remains
  unavailable.
- The effective-stress path is calibrated against three biased DSS histories;
  the density interpolation still needs untouched biased tests at other CSR,
  density, and confinement levels.
- Unbiased and loose histories retain their production limitations by design.
- No native OpenSees/Hercules port, implicit tangent, complete 3-D mesh test,
  or runtime qualification has been performed.

## Reproducibility files

- `RIVASandPhaseTransformation.py`: research kernel.
- `RIVASandPhaseTransformationTest.py`: disabled equivalence, restart,
  compatibility, finite-state, and substep check.
- `RIVASandPhaseTransformationAudit.py`: strict eight-case adapter.
- `plot_phase_transformation_comparison.py`: comparison figures.
- `results/phase_transformation_density_extended_final_full`: 32-by-4 histories
  and metrics.
- `results/phase_transformation_density_extended_final_64x2`: objectivity
  histories and metrics.
- `results/phase_transformation_density_extended_comparison`: final PNG/PDF
  figures, including isolated cycle-11 effective-stress paths, plus the
  machine-readable `effective_path_metrics.csv`.
