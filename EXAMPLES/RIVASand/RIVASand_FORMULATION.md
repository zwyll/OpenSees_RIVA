# RIVA-Sand: formulation and OpenSees material specification

RIVA expands to **R**eversible--**I**rreversible **V**olumetric
**A**ccumulation. `RIVASand` is the identifier used in source code and Tcl.

## 1. Scope and frozen status

This document specifies the equations embedded in the OpenSees
`RIVASand` material. The kernel is a
single, inheritance-free expansion of the calibrated V1--V8 Python model.  The
112 constitutive parameter defaults are identical to `pj_model_v8.reference`.
Four additional quantities control event detection at the host-step boundary;
they are numerical controls, not calibration parameters.

The V8 constitutive equations are frozen. The OpenSees `RIVASand`
interface only replaces ten values already present in the V8 parameter record;
it adds no equation and is not the V8.1 `h_p0` experiment. The only deliberate
algorithmic port change is the reversal detector in Section 12. Setting
`objective_reversal_enabled=False` restores the original substep-level detector
and reproduces the frozen implementation to roundoff.

V8 is an effective-stress, strain-driven material model.  In the DSS surrogate,
undrained loading is imposed with zero total volumetric strain and pore-pressure
buildup is represented by loss of effective mean pressure through volumetric
compatibility.  Do not combine this internal pressure-generation mechanism with
a second pore-pressure law without deriving a consistent coupled formulation.

## 2. Conventions and invariants

Stress and strain are tension-positive.  Effective compression therefore has
negative normal stress.  For effective Cauchy stress `sigma`:

```text
p       = -tr(sigma)/3                         (compression-positive)
s       = sigma + p I                          (deviatoric stress)
q       = sqrt((3/2) s:s)
alpha   = s/p                                  (stress-ratio tensor)
||A||   = sqrt(A:A)
```

All equations use true tensor shear components. OpenSees supplies engineering
shear strain `gamma_ij`; the adapter converts it at the interface using
`epsilon_ij = gamma_ij/2` for `i != j`.

## 3. Pressure-dependent elasticity and state surfaces

Let `p+ = max(p, p_min)`.  The reference elasticity and bounding-surface
hardening are

```text
E(p) = E_ref (p+/p_ref)^n_G
G(p) = E(p) / [2(1+nu)]
K(p) = E(p) / [3(1-2nu)]
h_p  = h (p+/p_ref)^(-q_H).
```

The pressure-dependent relative state is

```text
Dr       = clip[(e_max-e)/(e_max-e_min), 0, 1]
c(p)     = clip[R / max(Q-ln(100 p+/p_atm), 0.10), 0, 1]
xi_r     = clip[Dr-c(p), -relative_state_limit, relative_state_limit]
M_b      = M exp(n_b xi_r)
M_d      = kd exp(-n_d xi_r)
alpha_b  = sqrt(2/3) M_b
eta_d    = sqrt(2/3) M_d.
```

The stress ratio is projected back to `||alpha|| <= alpha_b` after each
explicit substep.

## 4. Bounding-surface mapping and plastic multiplier

The image point lies on the spherical bounding surface:

```text
a       = alpha - alpha0
alpha*  = alpha + beta a
||alpha*||^2 = (2/3) M_b^2
n       = alpha*/||alpha*||.
```

Thus `beta` is the positive root of

```text
(a:a) beta^2 + 2(a:alpha) beta
                    + alpha:alpha - (2/3)M_b^2 = 0.
```

If the quadratic is degenerate, its discriminant is negative, or its result is
not finite, `beta=beta0`.  The accepted value is bounded below by `1e-6`.

For a strain increment `Delta epsilon`, define

```text
Delta epsilon_v = tr(Delta epsilon)
Delta e         = Delta epsilon - (Delta epsilon_v/3) I
H                = p h_eff beta^m
den              = 2G + (2/3)H - K_c D (alpha:n)
num              = 2G (Delta e:n) + K_c Delta epsilon_v (alpha:n)
Delta lambda     = max(0, num/max[den, 2G denominator_floor_ratio]).
```

`K_c=K` normally and `K_c=0` when the compatibility pressure floor is already
active.  A floor use increments `denominator_floor_hits`.  The deviatoric trial
update is

```text
s_new = s_old + 2G(Delta e - Delta lambda n).
```

`h_eff` and the effective `G` include the cyclic, bias, and density corrections
defined in Sections 8--10.

## 5. Reversal overshoot rule

At an accepted reversal, with accumulated equivalent plastic increment
`ep_eq_since_reversal`:

```text
r       = ep_eq_since_reversal / overshoot_strain
omega   = 0                    if r > 1
          1-r                  otherwise
alpha0  = alpha_old            if ||alpha0_old|| = 0
          omega alpha01_old + (1-omega) alpha_old otherwise
alpha01 = alpha0_old
beta    = beta0
ep_eq_since_reversal = 0.
```

The physical event trigger used for the port is given in Section 12.

## 6. Reversible/irreversible volumetric compatibility

The dilatancy projection is `eta=alpha:n`.  The irreversible part is

```text
c_d = max(eta_d-eta, 0)
A_z = zeta [1 + max(z:n, 0)]
r_b = [sqrt(2/3)M_b-eta] / max(beta, 1e-12)
c_in = C_in_ratio sqrt(2/3)M_b
g_c = (r_b+c_in)^2 / (r_b^2+C_D)        when the gate is enabled
D_ir = A_z g_c c_d exp[-irreversible_decay max(-eps_v_ir,0)].
```

After calculation, `D_ir` is multiplied by the amplitude/state/bias contraction
factor in Sections 7 and 9.  The reversible part is

```text
D_re = -zeta reversible_ratio (eta-eta_d)       if eta > eta_d
D_re = min(reversible_release_rate eps_v_re,
           max(eta_d,1e-12))                    if eta <= eta_d and eps_v_re>0
D_re = 0                                        otherwise.
```

The compatible volumetric state advances as

```text
eps_v_ir,new = eps_v_ir,old - Delta lambda max(D_ir,0)
eps_v_re,new = max(eps_v_re,old - Delta lambda D_re, 0)
eps_v_total,new = eps_v_total,old + Delta epsilon_v
eps_v_c,new = eps_v_total,new - eps_v_ir,new - eps_v_re,new.
```

The invariant checked by `compatibility_residual` is

```text
eps_v_total - eps_v_c - eps_v_ir - eps_v_re = 0.
```

For `n_G != 1`, pressure and confining strain are mapped exactly through the
pressure-dependent bulk modulus:

```text
eps_v_c(p) = p_ref^n_G/[K_ref(1-n_G)]
             [p_anchor^(1-n_G)-p^(1-n_G)]

p(eps_v_c) = {p_anchor^(1-n_G)
              -(1-n_G)K_ref eps_v_c/p_ref^n_G}^{1/(1-n_G)}.
```

For `n_G=1`, use

```text
eps_v_c(p) = (p_ref/K_ref) ln(p_anchor/p)
p(eps_v_c) = p_anchor exp(-K_ref eps_v_c/p_ref).
```

Pressure is limited to `p_min`; `pressure_floor_hits` records each activation.

## 7. Cyclic amplitude and density-dependent contraction knee

At a reversal, let `s_r` be the deviatoric stress at the start of the accepted
reversal increment.  With the preceding reversal stress `s_last`:

```text
A = ||s_r-s_last|| / (d p_anchor),
d = 1 for the first reversal and 2 thereafter.
```

The density-adjusted knee uses the initial relative state `xi_r0`:

```text
Delta xi = xi_r0-reference_relative_state
k_eff = clip[low_amplitude_knee_ratio exp(k Delta xi),
             effective_knee_minimum, 1]
k = dense_knee_state_exponent for Delta xi >= 0,
    loose_knee_state_exponent otherwise.
```

The amplitude multiplier is

```text
x   = max(A/cyclic_amplitude_reference, 1e-12)
F_A = clip[x^cyclic_amplitude_exponent
           min(1,x/k_eff)^low_amplitude_exponent,
           amplitude_factor_minimum, amplitude_factor_maximum].
```

The state contraction factor is

```text
F_s = clip[exp(-b Delta xi), state_factor_minimum, state_factor_maximum]
b   = dense_state_contraction_exponent for Delta xi >= 0,
      loose_state_contraction_exponent otherwise.
```

The irreversible contraction multiplier is

```text
F_con = F_A F_s [1+bias_contraction_scale B^bias_contraction_exponent],
```

where projected bias `B` is defined below.

## 8. Cyclic modulus and base hardening corrections

After the required reversal count:

```text
x_p = clip(p/p_anchor,0,1)^cyclic_flow_pressure_exponent
F_G,cyc = 1-cyclic_shear_modulus_reduction x_p
F_H,cyc = 1+cyclic_hardening_boost x_p.
```

Before that count, both factors are one.

## 9. Static bias and ratcheting

At dynamic activation, V8 stores the normalized static-bias tensor

```text
b = (s_equilibrated-s_zero_bias_reference)/p_anchor.
```

Once a cyclic excursion direction `c` exists, the projected bias is

```text
B = |b:c|/||c||.
```

The bias hardening boost is

```text
A_b = bias_amplitude_ratio A
margin   = max(B-A_b,0)
crossing = max(A_b/B-1,0)
H_b = bias_hardening_intercept B^bias_intercept_exponent
      exp(-bias_crossing_decay crossing)
      + bias_hardening_scale margin^bias_margin_exponent.
```

It modifies hardening by

```text
F_H,bias = 1 + H_b C_p clip(p/p_anchor,0,1)^bias_pressure_exponent
C_p = clip(bias_reference_pressure/p_anchor,0.25,4)
      ^bias_confinement_exponent.
```

For the shear-ratchet branch, let `S(x)=clip(x,0,1)^2[3-2clip(x,0,1)]`:

```text
a_A = S[(A-A_on)/(A_full-A_on)]
a_R = 1-S[(A/B-r_full)/(r_cutoff-r_full)]
activity = a_A a_R
capacity = bias_ratchet_limit activity
           (bias_ratchet_reference_bias/B)^bias_ratchet_bias_exponent
rate = bias_ratchet_rate
       max(p_anchor/bias_reference_pressure-1,0)^bias_ratchet_pressure_exponent
       (bias_ratchet_reference_bias/B)^bias_ratchet_bias_exponent
       max(1-ratchet_strain/capacity,0)
Delta r = min(rate Delta lambda, capacity-ratchet_strain).
```

When active, the backbone update is repeated with
`Delta epsilon_eff = Delta epsilon-Delta r c_sign`, where `c_sign` is the unit
cyclic direction aligned with the static bias.

## 10. V8 density shakedown: decoupled hardening and compliance

Recover the initial relative density and the calibration-anchor density:

```text
Dr0 = clip[xi_r0 + R/max(Q-ln(100 p_anchor/p_atm),0.10),0,1]
Dr_ref = (e_max-state_shakedown_reference_void_ratio)/(e_max-e_min)
d = Dr0-Dr_ref.
```

After activation and the minimum reversal count, the raw V7 factor is

```text
F_raw = clip[1 + state_shakedown_scale
             {1-exp[-(d/w)^2]}
             exp(state_shakedown_state_sensitivity d)
             (B/state_shakedown_bias_reference)^state_shakedown_bias_exponent
             clip(p/p_anchor,0,1)^state_shakedown_pressure_exponent,
             1, state_shakedown_factor_maximum],
```

where `w=state_shakedown_anchor_width`.  It is one when disabled, before the
minimum reversal count, or at zero projected bias.

For the dense branch:

```text
w_d = 1-exp[-(max(d,0)/w)^2]
release = state_shakedown_dense_hardening_scale
          exp[-state_shakedown_dense_hardening_amplitude_decay
              max(A-state_shakedown_dense_hardening_amplitude_onset,0)]
F_H,state = clip[1+(F_raw-1){1+w_d(release-1)}, 1, F_max].
```

Compliance is attenuated independently:

```text
F_comp = 1+(F_raw-1)
         exp[-state_shakedown_dense_compliance_decay max(d,0)^2]
F_G,state = clip[F_comp^(-state_shakedown_compliance_exponent),
                 state_shakedown_shear_multiplier_minimum,1].
```

The final moduli and hardening used in Section 4 are

```text
G_eff = G F_G,cyc F_G,state
K_eff = K
h_eff = h_p F_H,cyc F_H,bias F_H,state.
```

## 11. Bias-dependent reversible pressure wave

After the first amplitude reversal, project the dynamic deviatoric stress on
the cyclic direction:

```text
s_static = s_geostatic + p_anchor b
phi = clip[(s-s_static):c_hat/(p_anchor A),-1,1].
```

Then

```text
v_osc = -bias_reversible_volume_amplitude
        (bias_reversible_volume_reference_bias/B)
          ^bias_reversible_volume_bias_exponent
        (p_anchor/bias_reference_pressure)
          ^bias_reversible_volume_pressure_exponent
        [1-exp(-N_r/bias_reversible_volume_buildup_reversals)] phi

v_mean = bias_reversible_mean_scale
         (bias_reversible_mean_transition_pressure-p_anchor)
         /bias_reference_pressure
         (B/bias_reversible_volume_reference_bias)
         [1-exp(-N_r/bias_reversible_mean_buildup_reversals)]

v_bias = v_osc+v_mean.
```

`N_r` is `amplitude_reversals`.  The pressure mapping uses
`eps_v_total = physical_eps_v_total+v_bias`, while void ratio advances only with
the physical volumetric strain.  Before each new mechanical update, the old
reversible bias volume is removed so it is not accumulated as permanent strain.

## 12. Host-increment reversal detector

The original detector evaluates

```text
(alpha_trial-alpha0):(alpha_trial-alpha_old) < 0
```

inside every constitutive substep.  Adaptive refinement or a change in the
number of fixed substeps can therefore create different event counts.  It is
available only as a regression mode.

The port detector instead executes once, before subdividing an OpenSees trial
increment:

```text
d_k = dev(Delta epsilon_k)/||dev(Delta epsilon_k)||
cos(theta) = d_k:d_(k-1)
reversal if cos(theta) <= reversal_direction_cosine
            and ||s_k-s_last|| >=
                reversal_stress_deadband_ratio p_anchor.
```

Increments with `||dev(Delta epsilon)|| <= reversal_strain_deadband` are ignored
and do not erase the previous direction.  At most one reversal is passed to the
first internal substep; all remaining internal substeps have event detection
disabled.  Detection starts only after `begin_cyclic_phase()`.

This makes event count independent of constitutive substepping and stable under
host timestep refinement when the strain turning points are resolved.  It does
not make an unresolved coarse history equivalent to a resolved history.

## 13. Complete explicit update order

For each OpenSees trial strain increment:

1. Validate a symmetric tensor increment and compute its host deviatoric
   direction.
2. Detect at most one host reversal using Section 12.
3. Split the increment into `nsub` equal fixed subincrements.
4. For each subincrement, remove the preceding reversible bias volume from the
   mechanical state.
5. Compute `p`, `s`, state-dependent `G_eff`, `K`, and `h_eff`.
6. Form the elastic trial state and apply the accepted reversal rule, if this is
   the first subincrement of a reversed host step.
7. Compute `Delta lambda`, deviatoric stress, compatible volumetric strains,
   pressure, void ratio, fabric, surface projection, `beta`, `n`, `D_ir`, and
   `D_re`.
8. At a reversal, update excursion amplitude, amplitude factor, last reversal
   stress, and reversal counters.
9. Evaluate ratchet increment from the trial; if nonzero, repeat the backbone
   update with the ratchet-corrected strain increment.
10. Update cyclic direction at a newly accepted reversal.
11. Compute reversible bias-volume target, rebuild compatible pressure, update
    physical void ratio, and scale `D_re` by the optional bias factor.
12. After the final substep, store the current nonzero host loading direction.

The method is forward Euler.  `adaptive`, `integration_rtol`,
`integration_atol`, `max_refinement_depth`, and `localize_events` remain in the
parameter record solely for exact frozen-parameter parity; this standalone
kernel deliberately ignores them.

## 14. State variables

| State | Type | Meaning / initialization |
|---|---:|---|
| `stress` | tensor | Effective stress `sigma`; supplied equilibrated value |
| `alpha` | tensor | `s/p`; computed from initial stress |
| `alpha0` | tensor | bounding-surface projection origin; zero |
| `alpha01` | tensor | previous projection origin used by overshoot; zero |
| `n` | tensor | mapped loading normal; zero |
| `fabric` | tensor | dilatancy fabric; zero |
| `D`, `D_ir`, `D_re` | scalar | total, irreversible, and reversible dilatancy; `D_ir=zeta eta_d`, `D_re=0` |
| `beta` | scalar | image-point mapping distance; `beta0` |
| `lambda_total` | scalar | cumulative plastic multiplier; zero |
| `ep_eq_since_reversal` | scalar | plastic multiplier since reversal; zero |
| `void_ratio` | scalar | physical void ratio; user input |
| `reversals` | integer | bounding-surface reversal count; zero |
| `pressure_floor_hits` | integer | pressure-floor diagnostic; zero |
| `denominator_floor_hits` | integer | plastic-denominator floor diagnostic; zero |
| `beta_fallbacks` | integer | quadratic mapping fallback count; zero |
| `eps_v_total` | scalar | compatible total volume, including reversible bias; zero |
| `eps_v_confining` | scalar | elastic/confining compatibility coordinate; zero |
| `eps_v_irreversible` | scalar | permanent contractive compatibility strain; zero |
| `eps_v_reversible` | scalar | recoverable dilative compatibility strain; zero |
| `pressure_anchor` | scalar | initial effective mean pressure; `p_initial` |
| `last_reversal_deviator` | tensor | deviatoric stress at last reversal; `s_initial` |
| `cyclic_amplitude` | scalar | latest normalized half/full excursion; zero |
| `amplitude_factor` | scalar | contraction amplitude multiplier; one |
| `amplitude_reversals` | integer | accepted excursion count; zero |
| `initial_relative_state` | scalar | frozen initial `xi_r`; computed |
| `state_contraction_factor` | scalar | frozen `F_s`; computed |
| `effective_knee_ratio` | scalar | frozen density-adjusted knee; computed |
| `geostatic_deviator` | tensor | deviatoric stress at initialization; `s_initial` |
| `static_bias_tensor` | tensor | normalized equilibrated bias; zero until activation |
| `cyclic_direction` | tensor | unit stress-excursion direction; zero |
| `static_bias_index` | scalar | `||s_equilibrated-s_reference||/p_anchor`; zero |
| `cyclic_phase_active` | logical | enables bias/reversal branches; false |
| `bias_ratchet_strain` | scalar | accumulated ratchet strain; zero |
| `physical_eps_v_total` | scalar | total volume excluding `v_bias`; zero |
| `bias_reversible_volume` | scalar | current reversible pressure-wave volume; zero |
| `last_host_deviatoric_strain_direction` | tensor | host reversal memory; zero |

All tensor states are symmetric.  A compact production layout may store each
as six components.  A clear recommended order is `[xx, yy, zz, xy, yz, xz]`;
the OpenSees adapter maps its native order to this order explicitly.

## 15. Parameters and frozen Ottawa F65 defaults

The units of `E_ref`, `p_ref`, `p_min`, `p_atm`, and pressure-valued bias
parameters must match the stress unit used by OpenSees. `h` is a dimensionless
hardening prefactor because the update multiplies it by effective pressure. The
frozen calibration below uses kPa for stress-like quantities.

### OpenSees material-property interface

The configurable command reads:

```text
nDMaterial RIVASand tag Dr M kd h m zeta eMax eMin Q R nG
```

`Dr` is dimensionless relative density and is converted at initialization by

```text
e0 = e_max - Dr (e_max-e_min).
```

The ten values following `Dr` are existing fields in the equations below:
`M` and `kd` set the base bounding and dilatancy stress ratios; `h` and `m`
control hardening; `zeta` scales irreversible and reversible dilatancy;
`e_max`, `e_min`, `Q`, and `R` define density/state mappings; and `n_G`
controls pressure dependence of elastic stiffness. `M` is supplied directly,
so the input does not also contain `phi_cs`. `q_H` remains independently
frozen at 0.35; it is not silently replaced by `1-n_G`.

OpenSees assigns one parameter row to each material tag; there is no automatic
Vs interpolation in this adapter. Inputs must be finite and satisfy
`0<=Dr<=1`, `M,kd,h,Q>0`,
`m,zeta,R>=0`, `e_max>e_min>=0`, and `0<=n_G<=1`. Because the otherwise
frozen V8 shakedown law retains `state_shakedown_reference_void_ratio=0.601`,
the custom interval must also satisfy `e_min<=0.601<=e_max`.

Every parameter not listed in the compact table—including `E_ref`, `nu`,
`q_H`, compatibility coefficients, amplitude scaling, static-bias response,
ratchet evolution, and reversal controls—retains the frozen default. A custom
row therefore needs calibration and sensitivity validation; the original
golden histories only certify the frozen/reference-valued row.

### Backbone, state, and legacy numerical record

| Parameter | Default | Role |
|---|---:|---|
| `E_ref` | 127339.75550887753 | reference Young modulus |
| `nu` | 0.30 | Poisson ratio |
| `h` | 122.44207260468994 | hardening prefactor |
| `m` | 0.945 | mapping-distance exponent |
| `M` | 1.25 | base bounding stress ratio |
| `kd` | 1.125 | base dilatancy stress ratio |
| `zeta` | 0.025 | dilatancy scale |
| `beta0` | 1000 | mapping fallback/reset |
| `overshoot_strain` | 5e-5 | reversal interpolation strain |
| `p_ref` | 101.3 | elastic pressure reference |
| `p_min` | 0.001 | effective-pressure floor |
| `n_G` | 0.65 | pressure exponent for stiffness |
| `q_H` | 0.35 | pressure exponent for hardening |
| `denominator_floor_ratio` | 0.02 | plastic denominator floor relative to `2G` |
| `state_enabled` | true | enable relative state |
| `e_max`, `e_min` | 0.78, 0.51 | limiting void ratios |
| `p_atm` | 101.3 | state pressure reference |
| `Q`, `R` | 10.0, 1.5 | critical-density relation |
| `n_b`, `n_d` | 0.05, 0.02 | bounding/dilatancy state exponents |
| `relative_state_limit` | 0.75 | absolute `xi_r` limit |
| `fabric_enabled` | true | enable fabric update |
| `c_z`, `z_max` | 600, 4 | fabric rate and limit |
| `contraction_gate_enabled` | true | enable inward gate |
| `C_D`, `C_in_ratio` | 0.10, 0.01 | contraction-gate controls |
| `adaptive` | true | retained; ignored in port kernel |
| `integration_rtol`, `integration_atol` | 0.003, 1e-8 | retained; ignored |
| `max_refinement_depth` | 6 | retained; ignored |
| `localize_events` | true | retained; ignored |

### Compatibility and cyclic amplitude

| Parameter | Default | Role |
|---|---:|---|
| `compatibility_enabled` | true | use volumetric pressure mapping |
| `reversible_enabled` | true | reversible dilatancy branch |
| `irreversible_enabled` | true | irreversible contraction branch |
| `reversible_ratio` | 7.10 | reversible dilation scale |
| `reversible_release_rate` | 300 | release rate |
| `irreversible_decay` | 0 | accumulated-contraction decay |
| `amplitude_scaling_enabled` | true | amplitude-dependent contraction |
| `cyclic_amplitude_reference` | 0.409 | amplitude normalization |
| `cyclic_amplitude_exponent` | 3.3 | main amplitude exponent |
| `low_amplitude_knee_ratio` | 0.85 | base knee ratio |
| `low_amplitude_exponent` | 6.0 | low-amplitude suppression |
| `amplitude_factor_minimum`, `amplitude_factor_maximum` | 0.05, 12.0 | amplitude-factor bounds |
| `state_contraction_enabled` | true | density-dependent contraction |
| `reference_relative_state` | 0.4405633323 | contraction/knee anchor |
| `dense_state_contraction_exponent` | 17.0 | dense-side exponent |
| `loose_state_contraction_exponent` | 6.5 | loose-side exponent |
| `state_factor_minimum`, `state_factor_maximum` | 0.005, 5.0 | state-factor bounds |
| `cyclic_flow_correction_enabled` | true | cyclic modulus/hardening branch |
| `cyclic_shear_modulus_reduction` | 0.85 | modulus reduction magnitude |
| `cyclic_hardening_boost` | 0 | cyclic hardening magnitude |
| `cyclic_flow_pressure_exponent` | 0.5 | pressure activity exponent |
| `cyclic_flow_minimum_reversals` | 1 | activation count |
| `state_dependent_knee_enabled` | true | density-adjusted knee |
| `loose_knee_state_exponent`, `dense_knee_state_exponent` | 3.7, 0 | knee state exponents |
| `effective_knee_minimum` | 0.20 | lower knee bound |

### Static-bias hardening and reversible pressure wave

| Parameter | Default | Role |
|---|---:|---|
| `static_bias_enabled` | true | enable all static-bias branches |
| `bias_hardening_intercept` | 50 | bias hardening intercept |
| `bias_intercept_exponent` | 2.3 | intercept bias exponent |
| `bias_crossing_decay` | 3.0 | decay after bias crossing |
| `bias_hardening_scale` | 335 | non-crossing margin scale |
| `bias_margin_exponent` | 1.5 | margin exponent |
| `bias_amplitude_ratio` | 1.0 | amplitude used in bias comparison |
| `bias_pressure_exponent` | 0.5 | current-pressure hardening exponent |
| `bias_reference_pressure` | 26.266666666666666 | bias pressure reference |
| `bias_confinement_exponent` | 0.85 | initial-confinement exponent |
| `bias_minimum_reversals` | 1 | bias branch activation count |
| `bias_contraction_scale`, `bias_contraction_exponent` | 0, 1 | optional contraction multiplier |
| `bias_reversible_scale`, `bias_reversible_exponent` | 0, 1 | optional `D_re` multiplier |
| `bias_reversible_volume_enabled` | true | pressure-wave volume branch |
| `bias_reversible_volume_amplitude` | 7.6e-5 | oscillatory volume amplitude |
| `bias_reversible_volume_reference_bias` | 0.646 | reference projected bias |
| `bias_reversible_volume_bias_exponent` | 2.0 | inverse-bias exponent |
| `bias_reversible_volume_pressure_exponent` | 0.22 | confinement exponent |
| `bias_reversible_volume_buildup_reversals` | 6.0 | oscillation buildup count |
| `bias_reversible_mean_scale` | 9e-5 | mean pressure-wave shift scale |
| `bias_reversible_mean_transition_pressure` | 40.0 | sign-change pressure |
| `bias_reversible_mean_buildup_reversals` | 12.0 | mean buildup count |

### Shear ratchet and V8 density shakedown

| Parameter | Default | Role |
|---|---:|---|
| `bias_ratchet_enabled` | true | enable biased shear ratchet |
| `bias_ratchet_rate` | 2.0 | ratchet rate per plastic multiplier |
| `bias_ratchet_limit` | 0.0215 | base strain capacity |
| `bias_ratchet_amplitude_onset`, `bias_ratchet_amplitude_full` | 0.30, 0.43 | amplitude activity window |
| `bias_ratchet_ratio_full`, `bias_ratchet_ratio_cutoff` | 0.70, 0.80 | amplitude/bias cutoff window |
| `bias_ratchet_reference_bias` | 0.646 | ratchet reference bias |
| `bias_ratchet_bias_exponent` | 0.62 | bias exponent |
| `bias_ratchet_pressure_exponent` | 1.0 | confinement exponent |
| `state_shakedown_enabled` | true | V7/V8 density branch |
| `state_shakedown_scale` | 5.0 | raw hardening magnitude |
| `state_shakedown_state_sensitivity` | 2.0 | density exponential |
| `state_shakedown_reference_void_ratio` | 0.601 | calibration anchor void ratio |
| `state_shakedown_anchor_width` | 0.02 | smooth anchor width |
| `state_shakedown_bias_reference` | 0.646 | shakedown reference bias |
| `state_shakedown_bias_exponent` | 0.5 | shakedown bias exponent |
| `state_shakedown_pressure_exponent` | 0 | pressure exponent |
| `state_shakedown_minimum_reversals` | 1 | activation count |
| `state_shakedown_factor_maximum` | 50 | hardening-factor ceiling |
| `state_shakedown_compliance_exponent` | 1.5 | compliance power |
| `state_shakedown_shear_multiplier_minimum` | 0.03 | shear-multiplier floor |
| `state_shakedown_dense_hardening_scale` | 1.0 | dense hardening scale |
| `state_shakedown_dense_hardening_amplitude_onset` | 0.82 | dense release onset |
| `state_shakedown_dense_hardening_amplitude_decay` | 3.5 | dense release rate |
| `state_shakedown_dense_compliance_decay` | 200 | dense compliance attenuation |

### Port-only event controls

| Parameter | Default | Role |
|---|---:|---|
| `objective_reversal_enabled` | true | host-step detector; false reproduces frozen legacy logic |
| `reversal_direction_cosine` | -0.20 | direction-change threshold |
| `reversal_strain_deadband` | 1e-12 | ignore nearly zero deviatoric increments |
| `reversal_stress_deadband_ratio` | 1e-4 | minimum stress excursion divided by `p_anchor` |

## 16. Initialization and OpenSees lifecycle

For ordinary element analysis:

1. Create the material in stage 0. It uses the frozen reference elastic
   tangent while the gravity/geostatic analysis establishes committed
   **effective** stress at every integration point.
2. Commit gravity, then call `updateMaterialStage -material tag -stage 1`.
   Each material copy initializes from its own committed stress. Its normal
   stress is the geostatic reference; committed shear becomes static bias.
3. During dynamics, OpenSees passes the total trial strain. The adapter forms
   the increment from the last committed strain, converts engineering shear,
   and calls the fixed-substep kernel. Start with `-nSub 1`; increase it only
   when a material-point substep study justifies the cost.
4. `commitState`, `revertToLastCommit`, `sendSelf`, and `recvSelf` preserve the
   complete state in Section 14.

The frozen constitutive guard is `p_min = 0.001 kPa` (scaled by
`-stressScale`). OpenSees also applies a separate tangent-only pressure floor
when evaluating the elastic tangent returned to the global solver. Its default
is `p_ref/200`; `-tangentPMin value` overrides it without changing the stress
update. `-pMin value` instead overrides the constitutive floor and therefore
must be treated as a material-model change in validation studies. Both option
values use the current OpenSees stress unit.

The stock historical `MaterialStageParameter` stopped after finding the first
matching element because older soil models stored stage globally. V8 stores
the activation state and geostatic anchor at each integration point, so stage
activation must register and update every matching element-held material copy.

For a direct material-point test, `-initialStress sxx syy szz sxy syz sxz`
starts directly in stage 1. Normal compression is negative.

If gravity is solved with another elastic/geostatic procedure, import its final
effective stress rather than replaying the gravity path through V8.  The helper
`initialize_equilibrated_state()` demonstrates this contract; it is not a
replacement for a spatial gravity solution.

In a coupled OpenSees `u-p` formulation, V8 owns effective skeleton stress and
its internal elastic/plastic volumetric compatibility. `BrickUP`, `BBarBrickUP`,
or `SSPbrickUP` independently owns the pressure degree of freedom through its
fluid mass balance, Darcy-flow, and displacement-pressure coupling terms.
Do not copy the standalone inferred `r_u` or V8 effective-pressure-loss
diagnostic into the BrickUP pressure state. Station `r_u` for a `u-p` run must
come from the solved excess pore pressure, `Delta p_f/sigma_v0'`.

## 17. Driver outputs and verification acceptance criteria

`RIVASand_material_point.tcl` imposes fixed host strain increments on a
homogeneous unit `bbarBrick` and reports:

```text
ru_vertical = 1 - (-sigma_zz)/sigma_v0'
ru_mean     = 1 - p/p0.
```

The focused tests require:

- all 112 frozen constitutive parameters to match exactly;
- all common frozen V8 state fields to match to `1e-12` when the legacy
  detector is selected;
- ten reversals for five resolved sinusoidal cycles with 1, 2, 4, or 8 internal
  substeps; and
- eight reversals for four cycles with 16, 32, 64, or 128 host points per cycle.

In addition, the linked OpenSees Tcl run is compared with
`golden_data/cyclic_zero_bias_reference.csv`. Its maximum accepted absolute
differences are `5e-6` in stress, `2e-7` in the skeleton pressure-loss
diagnostic, and `1e-12` in the compatibility residual; reversal counts must be
identical. The small stress tolerance accommodates the brick equilibrium solve
while remaining several orders below the constitutive response scale.

These checks establish implementation equivalence and event-count objectivity;
they do not replace the existing PRJ-3484/PRJ-4666 validation.  V8's known dense
biased-test limitation remains part of the frozen calibration and must not be
silently retuned during the OpenSees translation.
