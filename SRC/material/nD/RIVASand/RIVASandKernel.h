/* -*- C++ -*- */
#ifndef OPENSEES_RIVA_SAND_KERNEL_H
#define OPENSEES_RIVA_SAND_KERNEL_H

/* Native, allocation-free RIVA-Sand constitutive kernel.
 *
 * OpenSees copy: the constitutive equations and operation order are identical
 * to the verified Hercules RIVA-Sand kernel (parameter SHA below). The OpenSees
 * NDMaterial adapter owns engineering-shear conversion, lifecycle, and
 * serialization; this header remains solver-independent.
 *
 * Source of truth: rivasand_port/model.py and RIVASAND_FORMULATION.md.
 * Tensor order is [xx, yy, zz, xy, yz, xz], with physical
 * (not engineering) shear components. Stress and strain are tension-positive.
 *
 * Coupled u-p element ownership is deliberately outside this allocation-free
 * kernel.  The pressure reconstructed from eps_v_confining is effective
 * skeleton mean stress p', including the response to RIVA-Sand plastic volumetric
 * strain; it is never written to the element fluid-pressure field.  In a
 * coupled element OpenSees independently solves fluid mass balance
 *
 *     d(zeta) = alpha*d(eps_v) + d(p_f)/M,
 *
 * and equilibrium uses sigma_total=sigma_effective-alpha*p_f*I.  Therefore
 * the standalone diagnostic 1-p'/p'_anchor must not be injected as p_f.
 */

#include <float.h>
#include <math.h>
#include <stdint.h>

typedef struct riva_tensor_t {
    double xx, yy, zz, xy, yz, xz;
} riva_tensor_t;
typedef struct riva_material_parameters_t {
    double h, m, M, kd, zeta, e_max, e_min, Q, R, n_G;
} riva_material_parameters_t;

#if defined(__CUDACC__)
#define RIVA_HD __host__ __device__
#else
#define RIVA_HD
#endif

#define RIVA_PARAMETER_COUNT 118
#define RIVA_CONSTITUTIVE_PARAMETER_COUNT 114
#define RIVA_STATE_VALUE_COUNT 95
#define RIVA_PARAMETER_SHA256 \
    "9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83"

typedef struct riva_parameters_t {
    double E_ref, nu, h, m, M, kd, zeta, beta0, overshoot_strain;
    double p_ref, p_min, p_residual, n_G, q_H, denominator_floor_ratio;
    int32_t geostatic_admission_enabled;
    int32_t state_enabled;
    double e_max, e_min, p_atm, Q, R, n_b, n_d, relative_state_limit;
    int32_t fabric_enabled;
    double c_z, z_max;
    int32_t contraction_gate_enabled;
    double C_D, C_in_ratio;
    int32_t adaptive;
    double integration_rtol, integration_atol;
    int32_t max_refinement_depth, localize_events;
    int32_t compatibility_enabled, reversible_enabled, irreversible_enabled;
    double reversible_ratio, reversible_release_rate, irreversible_decay;
    int32_t amplitude_scaling_enabled;
    double cyclic_amplitude_reference, cyclic_amplitude_exponent;
    double low_amplitude_knee_ratio, low_amplitude_exponent;
    double amplitude_factor_minimum, amplitude_factor_maximum;
    int32_t state_contraction_enabled;
    double reference_relative_state, dense_state_contraction_exponent;
    double loose_state_contraction_exponent, state_factor_minimum;
    double state_factor_maximum;
    int32_t cyclic_flow_correction_enabled;
    double cyclic_shear_modulus_reduction, cyclic_hardening_boost;
    double cyclic_flow_pressure_exponent;
    int32_t cyclic_flow_minimum_reversals;
    int32_t state_dependent_knee_enabled;
    double loose_knee_state_exponent, dense_knee_state_exponent;
    double effective_knee_minimum;
    int32_t static_bias_enabled;
    /* Research diagnostics (attribution of the sub-triggering damping hole).
     * Both default 0 = frozen behavior bit-for-bit. diag_no_bias_hardening
     * zeroes ONLY the bias-hardening boost inside riva_hardening_for_state;
     * diag_no_cyclic_flow makes riva_cyclic_flow return unity factors. */
    int32_t diag_no_bias_hardening, diag_no_cyclic_flow;
    /* L3 research: continuous plastic-activity scaling of the two
     * sub-triggering suppressors (attribution 2026-08-25: the bias-hardening
     * boost supplies the "elastic" half of the zero-dissipation ribbons, the
     * cyclic-flow shear cut the "soft" half; zeroing both overshoots to
     * xi ~ 28%). activity = min(1, ep_half_last / l3_ep_ref), latched per
     * half-cycle at reversal registration, so there is no per-step feedback.
     * Each suppressor scales as floor + (1-floor)*activity: fully calibrated
     * behavior at high plastic activity, partially released ribbons at low.
     * l3_ep_ref <= 0 disables the mechanism bit-for-bit (frozen). */
    double l3_ep_ref, l3_boost_floor, l3_cut_floor;
    /* V11 research (CycLiq-style memory-surface reload plasticity): when
     * enabled, (a) the cyclic-flow shear cut is off, (b) the bias-hardening
     * boost is REPLACED by a stress-ratio-memory factor
     * (1/clip(eta_m/mb, v11_eta_floor, 1))^v11_k: virgin/biased states are
     * protected like the boost protected them, but the factor DECAYS to 1
     * once strong shaking has driven eta near the bound, so reload branches
     * regain plasticity (hysteretic damping + Dre exchange). Default
     * disabled = frozen bit-for-bit. */
    int32_t v11_enabled;
    double v11_k, v11_eta_floor;
    /* V11 A2: floor on the amplitude factor AS APPLIED TO D_ir only.
     * The 0.05 amplitude floor starves plastic contraction at small cyclic
     * amplitudes (field CSR ~0.1), which is what the cyclic-flow shear cut
     * was compensating for as a strain-amplifying pump. CycLiqCPSP's Dir is
     * amplitude-agnostic; v11_dir_floor = 1.0 reproduces that, smaller
     * values keep part of the calibrated amplitude shape. */
    double v11_dir_floor;
    /* V11 A4 (CycLiqCPSP Dir transplant, the two decay terms that shape the
     * N-dependence of biased triggering):
     * - v11_rdr: within-branch "gammamonos" decay, w = (rdr/(rdr+ep_eq))^2
     *   with ep_eq = ep_eq_since_reversal (resets at each registered
     *   reversal). Fresh contraction burst per branch, fading with monotonic
     *   plastic strain - CycLiq's rdr term.
     * - v11_eta_ir: saturation brake exp(-eta_ir*|eps_v_irreversible|) -
     *   CycLiq's exp(-eta*epsvir): accumulated irreversible contraction
     *   progressively shuts further contraction down (the high-N brake).
     * Both default 0 = off (bit-exact even under v11). */
    double v11_rdr, v11_eta_ir;
    /* V12: wholesale replacement of the amplitude overlay on D_ir. The
     * frozen dilatancy core ALREADY carries the CycLiq structure - the
     * (rd - r):n projection (distance = sqrt(2/3) md - alpha:n), fabric
     * amplification, the C_D contraction gate, and an eps_v_ir decay. The
     * amplitude machinery is only a multiplicative overlay
     * (riva_irreversible_factor). Under v12 that overlay becomes a pure
     * state-coupled chain: state_contraction * v12_dir_scale * exp(chi)
     * with chi = min(1, v12_chi * eps_v_re/|eps_v_ir|) (CycLiqCPSP's
     * reversible/irreversible exchange coupling), plus the shared A4 knobs
     * (v11_eta_ir saturation, v11_rdr gammamonos). NO amplitude factor, NO
     * bias-contraction boost: amplitude-gradedness and K_alpha come from
     * the projection geometry itself. Pair with -v11 for the hardening
     * side. Default off = frozen bit-for-bit. */
    int32_t v12_enabled;
    double v12_dir_scale, v12_chi;
    double bias_hardening_intercept, bias_intercept_exponent;
    double bias_crossing_decay, bias_hardening_scale, bias_margin_exponent;
    double bias_amplitude_ratio, bias_pressure_exponent;
    double bias_reference_pressure, bias_confinement_exponent;
    int32_t bias_minimum_reversals;
    double bias_contraction_scale, bias_contraction_exponent;
    double bias_reversible_scale, bias_reversible_exponent;
    int32_t bias_reversible_volume_enabled;
    double bias_reversible_volume_amplitude;
    double bias_reversible_volume_reference_bias;
    double bias_reversible_volume_bias_exponent;
    double bias_reversible_volume_pressure_exponent;
    double bias_reversible_volume_buildup_reversals;
    /* V9 PROPORTIONAL amplitude consistency (breaks V8 bit-lock ??oracle
       re-frozen): the reversible-volume response scales LINEARLY with the
       plastic multiplier of the last completed half-cycle,
       min(1, ep_half_last/ep_cal) ??PM4Sand/SANISAND-style: elastic or
       weakly-plastic cycles drive proportionally little volume. ep_cal sits
       below the normal calibrated half-cycles (1.6e-3..6.7e-3), so only
       weak-plasticity half-cycles are rescaled. */
    double bias_reversible_volume_ep_cal;
    /* V9 bias-factor knee: (ref/bias)^2 diverges as bias -> 0, overdriving
       gentle-slope states ~11x above the calibrated envelope. At and above
       the knee the frozen expression applies; below it the factor fades
       linearly to zero (peak (ref/knee)^2 = 4.64 at knee 0.30, the same
       scale as the calibrated bias025 case). */
    double bias_reversible_volume_bias_knee;
    /* Low-pressure fade. The calibrated envelope anchors at ~26 kPa; at low
       pressure anchors the fixed volumetric amplitude becomes comparable to
       the confinement (small bulk modulus), so the whole response fades with
       smoothstep(p_anchor/pressure_fade). Exactly 1x the frozen expression
       at calibrated anchors (>= pressure_fade). */
    double bias_reversible_volume_pressure_fade;
    double bias_reversible_mean_scale;
    double bias_reversible_mean_transition_pressure;
    double bias_reversible_mean_buildup_reversals;
    int32_t bias_ratchet_enabled;
    double bias_ratchet_rate, bias_ratchet_limit;
    double bias_ratchet_amplitude_onset, bias_ratchet_amplitude_full;
    double bias_ratchet_ratio_full, bias_ratchet_ratio_cutoff;
    double bias_ratchet_reference_bias, bias_ratchet_bias_exponent;
    double bias_ratchet_pressure_exponent;
    int32_t state_shakedown_enabled;
    double state_shakedown_scale, state_shakedown_state_sensitivity;
    double state_shakedown_reference_void_ratio;
    double state_shakedown_anchor_width, state_shakedown_bias_reference;
    double state_shakedown_bias_exponent, state_shakedown_pressure_exponent;
    int32_t state_shakedown_minimum_reversals;
    double state_shakedown_factor_maximum;
    double state_shakedown_compliance_exponent;
    double state_shakedown_shear_multiplier_minimum;
    double state_shakedown_dense_hardening_scale;
    double state_shakedown_dense_hardening_amplitude_onset;
    double state_shakedown_dense_hardening_amplitude_decay;
    double state_shakedown_dense_compliance_decay;
    int32_t objective_reversal_enabled;
    double reversal_direction_cosine, reversal_strain_deadband;
    double reversal_stress_deadband_ratio;
    /* Admit a committed stress ratio inherited from an external (e.g. elastic
     * geostatic) phase that lies beyond the bounding surface: the projection
     * limit becomes max(bound, |alpha_committed|), so the inherited overshoot
     * is never snap-projected away in a single call (a finite stress jump for
     * an infinitesimal strain increment, which breaks implicit equilibrium),
     * while any FURTHER excursion beyond the inherited ratio is still capped.
     * The overshoot decays only through the constitutive response itself.
     * Default 0 preserves the frozen behavior bit-for-bit. */
    int32_t admit_inherited_overbound;
    /* Opt-in regularization: floor on the bounding-distance variable beta as
     * used in the plastic hardening modulus (hardening = p*h_eff*beta^m).
     * States at/beyond the bounding surface otherwise carry beta ~ 0 and a
     * near-zero plastic modulus, which makes implicit equilibrium singular
     * (SSI foundation-edge states under sustained static bias). The floor
     * applies ONLY at the hardening use-site; beta's stored evolution is
     * untouched. Scale: beta0 = 1e3 at each reversal; a floor of ~1 gives
     * percent-level residual stiffness. Default 0 preserves the frozen
     * behavior bit-for-bit. */
    double beta_floor;
    /* Opt-in PM4Sand-style additive hardening reserve: hardening becomes
     * p*h_eff*(beta^m + beta_reserve^m), but ONLY at material points that
     * have not yet registered a reversal (reversals == 0) ??i.e. the
     * monotonic stage-in segment where SSI foundation states approach the
     * bounding surface and the plastic modulus would otherwise vanish
     * (a hard floor fixes the value but keeps a kink; the additive form is
     * smooth). After the first reversal the exact frozen formula applies,
     * so cyclic calibration is untouched by construction. Default 0
     * preserves the frozen behavior bit-for-bit. */
    double beta_reserve;
    /* Opt-in ceiling on beta AT THE HARDENING USE-SITE ONLY (research:
     * sub-triggering hysteretic damping). Right after a reversal the mapping
     * distance beta is huge (beta0 = 1e3 at re-registration), so the plastic
     * modulus p*h_eff*beta^m is effectively infinite and every unload-reload
     * branch is elastic: zero-area loops, xi ~ 0.5% at moderate strain where
     * element tests show ~9%. Capping beta bounds the early-branch plastic
     * modulus so each branch dissipates a controlled amount, while the
     * near-bound region (beta << cap) that controls strength and triggering
     * is untouched. beta's stored evolution, the mapping, and dilatancy all
     * see the true beta. Default 0 preserves the frozen behavior
     * bit-for-bit. */
    double beta_cap;
    /* Opt-in ceiling on mean effective stress. Garbage Newton trial strains
     * (equilibrium exploration, not physics) otherwise produce astronomic
     * pressures that push the linear solver into numerically fragile
     * territory. Set far above any physical state (e.g. 50x overburden) so
     * it binds only on doomed trials. Default 0 disables (frozen). */
    double p_max;
} riva_parameters_t;

typedef struct riva_state_t {
    riva_tensor_t stress, alpha, alpha0, alpha01, n, fabric;
    double D, beta, lambda_total, ep_eq_since_reversal, void_ratio;
    int64_t reversals, pressure_floor_hits, denominator_floor_hits;
    int64_t beta_fallbacks;
    double eps_v_total, eps_v_confining, eps_v_irreversible;
    double eps_v_reversible, pressure_anchor, D_ir, D_re;
    riva_tensor_t last_reversal_deviator;
    double cyclic_amplitude, amplitude_factor;
    int64_t amplitude_reversals;
    double initial_relative_state, state_contraction_factor;
    double effective_knee_ratio;
    riva_tensor_t geostatic_deviator, static_bias_tensor, cyclic_direction;
    double static_bias_index;
    int32_t cyclic_phase_active;
    double bias_ratchet_strain, physical_eps_v_total;
    double bias_reversible_volume;
    riva_tensor_t last_host_deviatoric_strain_direction;
    /* plastic multiplier accumulated over the last COMPLETED half-cycle
       (captured from ep_eq_since_reversal at reversal registration) */
    double ep_half_last;
    /* V11 research: memory of the largest normalized stress ratio
     * eta = sqrt(3/2)|alpha| ever committed (CycLiqCPSP's eta_m). Drives
     * the memory-surface hardening factor that replaces the bias boost
     * under -v11. Tracked unconditionally (cheap), CONSUMED only when
     * v11_enabled - the frozen path never reads it. */
    double eta_m;
    int32_t initialized;
} riva_state_t;

typedef struct riva_update_info_t {
    int32_t accepted_substeps;
    int32_t reversal_registered;
} riva_update_info_t;

RIVA_HD static inline double riva_min(double a, double b) { return a < b ? a : b; }
RIVA_HD static inline double riva_max(double a, double b) { return a > b ? a : b; }
RIVA_HD static inline double riva_clip(double x, double a, double b)
{ return riva_min(riva_max(x, a), b); }
RIVA_HD static inline riva_tensor_t riva_zero(void)
{ riva_tensor_t a = {0,0,0,0,0,0}; return a; }
RIVA_HD static inline riva_tensor_t riva_add(riva_tensor_t a, riva_tensor_t b)
{ riva_tensor_t c={a.xx+b.xx,a.yy+b.yy,a.zz+b.zz,a.xy+b.xy,a.yz+b.yz,a.xz+b.xz}; return c; }
RIVA_HD static inline riva_tensor_t riva_sub(riva_tensor_t a, riva_tensor_t b)
{ riva_tensor_t c={a.xx-b.xx,a.yy-b.yy,a.zz-b.zz,a.xy-b.xy,a.yz-b.yz,a.xz-b.xz}; return c; }
RIVA_HD static inline riva_tensor_t riva_scale(riva_tensor_t a, double x)
{ riva_tensor_t c={a.xx*x,a.yy*x,a.zz*x,a.xy*x,a.yz*x,a.xz*x}; return c; }
RIVA_HD static inline riva_tensor_t riva_iso(double x)
{ riva_tensor_t a={x,x,x,0,0,0}; return a; }
RIVA_HD static inline double riva_trace(riva_tensor_t a) { return a.xx+a.yy+a.zz; }
RIVA_HD static inline double riva_ddot(riva_tensor_t a, riva_tensor_t b)
{ return a.xx*b.xx+a.yy*b.yy+a.zz*b.zz+2.0*(a.xy*b.xy+a.yz*b.yz+a.xz*b.xz); }
RIVA_HD static inline double riva_norm(riva_tensor_t a)
{ return sqrt(riva_max(riva_ddot(a,a),0.0)); }
RIVA_HD static inline riva_tensor_t riva_dev(riva_tensor_t a)
{ return riva_sub(a,riva_iso(riva_trace(a)/3.0)); }
RIVA_HD static inline double riva_pressure(riva_tensor_t stress)
{ return -riva_trace(stress)/3.0; }
/* Research low-confinement translation.  With p_residual=0 this is exactly
 * the frozen RIVA-Sand pressure.  A positive value translates the cone apex
 * without changing the physical effective-stress tensor returned to the
 * host. */
RIVA_HD static inline double riva_cone_pressure(
    const riva_parameters_t *p,riva_tensor_t stress)
{ return riva_pressure(stress)+p->p_residual; }
RIVA_HD static inline double riva_physical_pressure(
    const riva_parameters_t *p,double cone_pressure)
{ return cone_pressure-p->p_residual; }
/* p_min is the minimum physical effective pressure.  The translated cone
 * therefore has its apex at p_residual+p_min, not at p_min.  This distinction
 * is immaterial when p_residual=0 and preserves the frozen RIVA-Sand oracle. */
RIVA_HD static inline double riva_cone_pressure_floor(
    const riva_parameters_t *p)
{ return p->p_residual+p->p_min; }
RIVA_HD static inline double riva_q(riva_tensor_t stress)
{ riva_tensor_t s=riva_dev(stress); return sqrt(riva_max(1.5*riva_ddot(s,s),0.0)); }
RIVA_HD static inline int riva_finite_tensor(riva_tensor_t a)
{ return isfinite(a.xx)&&isfinite(a.yy)&&isfinite(a.zz)&&isfinite(a.xy)&&isfinite(a.yz)&&isfinite(a.xz); }

/* stress_scale converts the frozen kPa calibration to the caller's stress
 * unit. Use 1 for kPa and 1000 for Pa. */
RIVA_HD static inline riva_parameters_t riva_reference_parameters(double stress_scale)
{
    riva_parameters_t p = {};
    p.E_ref=127339.75550887753*stress_scale; p.nu=0.30;
    p.h=122.44207260468994; p.m=0.945; p.M=1.25; p.kd=1.125;
    p.zeta=0.025; p.beta0=1.0e3; p.overshoot_strain=5.0e-5;
    p.p_ref=101.3*stress_scale; p.p_min=1.0e-3*stress_scale;
    p.p_residual=0.0; p.geostatic_admission_enabled=0;
    p.n_G=0.65; p.q_H=0.35; p.denominator_floor_ratio=0.02;
    p.state_enabled=1; p.e_max=0.78; p.e_min=0.51;
    p.p_atm=101.3*stress_scale; p.Q=10.0; p.R=1.5;
    p.n_b=0.05; p.n_d=0.02; p.relative_state_limit=0.75;
    p.fabric_enabled=1; p.c_z=600.0; p.z_max=4.0;
    p.contraction_gate_enabled=1; p.C_D=0.10; p.C_in_ratio=0.01;
    p.adaptive=1; p.integration_rtol=3.0e-3; p.integration_atol=1.0e-8;
    p.max_refinement_depth=6; p.localize_events=1;
    p.compatibility_enabled=1; p.reversible_enabled=1;
    p.irreversible_enabled=1; p.reversible_ratio=7.10;
    p.reversible_release_rate=300.0; p.irreversible_decay=0.0;
    p.amplitude_scaling_enabled=1; p.cyclic_amplitude_reference=0.409;
    p.cyclic_amplitude_exponent=3.3; p.low_amplitude_knee_ratio=0.85;
    p.low_amplitude_exponent=6.0; p.amplitude_factor_minimum=0.05;
    p.amplitude_factor_maximum=12.0; p.state_contraction_enabled=1;
    p.reference_relative_state=0.4405633323;
    p.dense_state_contraction_exponent=17.0;
    p.loose_state_contraction_exponent=6.5; p.state_factor_minimum=0.005;
    p.state_factor_maximum=5.0; p.cyclic_flow_correction_enabled=1;
    p.cyclic_shear_modulus_reduction=0.85; p.cyclic_hardening_boost=0.0;
    p.cyclic_flow_pressure_exponent=0.5; p.cyclic_flow_minimum_reversals=1;
    p.state_dependent_knee_enabled=1; p.loose_knee_state_exponent=3.7;
    p.dense_knee_state_exponent=0.0; p.effective_knee_minimum=0.20;
    p.static_bias_enabled=1;
    p.diag_no_bias_hardening=0; p.diag_no_cyclic_flow=0;
    p.l3_ep_ref=0.0; p.l3_boost_floor=0.0; p.l3_cut_floor=0.0;
    p.v11_enabled=0; p.v11_k=2.0; p.v11_eta_floor=0.35;
    p.v11_dir_floor=1.0; p.v11_rdr=0.0; p.v11_eta_ir=0.0;
    p.v12_enabled=0; p.v12_dir_scale=1.0; p.v12_chi=0.0;
    p.bias_hardening_intercept=50.0;
    p.bias_intercept_exponent=2.3; p.bias_crossing_decay=3.0;
    p.bias_hardening_scale=335.0; p.bias_margin_exponent=1.5;
    p.bias_amplitude_ratio=1.0; p.bias_pressure_exponent=0.5;
    p.bias_reference_pressure=26.266666666666666*stress_scale;
    p.bias_confinement_exponent=0.85; p.bias_minimum_reversals=1;
    p.bias_contraction_scale=0.0; p.bias_contraction_exponent=1.0;
    p.bias_reversible_scale=0.0; p.bias_reversible_exponent=1.0;
    p.bias_reversible_volume_enabled=1;
    p.bias_reversible_volume_amplitude=7.6e-5;
    p.bias_reversible_volume_reference_bias=0.646;
    p.bias_reversible_volume_bias_exponent=2.0;
    p.bias_reversible_volume_pressure_exponent=0.22;
    p.bias_reversible_volume_buildup_reversals=6.0;
    p.bias_reversible_volume_ep_cal=1.0e-3;
    p.bias_reversible_volume_bias_knee=0.30;
    p.bias_reversible_volume_pressure_fade=12.0*stress_scale;
    p.bias_reversible_mean_scale=9.0e-5;
    p.bias_reversible_mean_transition_pressure=40.0*stress_scale;
    p.bias_reversible_mean_buildup_reversals=12.0;
    p.bias_ratchet_enabled=1; p.bias_ratchet_rate=2.0;
    p.bias_ratchet_limit=0.0215; p.bias_ratchet_amplitude_onset=0.30;
    p.bias_ratchet_amplitude_full=0.43; p.bias_ratchet_ratio_full=0.70;
    p.bias_ratchet_ratio_cutoff=0.80; p.bias_ratchet_reference_bias=0.646;
    p.bias_ratchet_bias_exponent=0.62; p.bias_ratchet_pressure_exponent=1.0;
    p.state_shakedown_enabled=1; p.state_shakedown_scale=5.0;
    p.state_shakedown_state_sensitivity=2.0;
    p.state_shakedown_reference_void_ratio=0.601;
    p.state_shakedown_anchor_width=0.02;
    p.state_shakedown_bias_reference=0.646;
    p.state_shakedown_bias_exponent=0.5;
    p.state_shakedown_pressure_exponent=0.0;
    p.state_shakedown_minimum_reversals=1;
    p.state_shakedown_factor_maximum=50.0;
    p.state_shakedown_compliance_exponent=1.5;
    p.state_shakedown_shear_multiplier_minimum=0.03;
    p.state_shakedown_dense_hardening_scale=1.0;
    p.state_shakedown_dense_hardening_amplitude_onset=0.82;
    p.state_shakedown_dense_hardening_amplitude_decay=3.5;
    p.state_shakedown_dense_compliance_decay=200.0;
    p.objective_reversal_enabled=1; p.reversal_direction_cosine=-0.20;
    p.reversal_strain_deadband=1.0e-12;
    p.reversal_stress_deadband_ratio=1.0e-4;
    p.admit_inherited_overbound=0;
    p.beta_floor=0.0;
    p.beta_reserve=0.0;
    p.beta_cap=0.0;
    p.p_max=0.0;
    return p;
}

RIVA_HD static inline riva_material_parameters_t riva_reference_material_parameters(
    const riva_parameters_t *p)
{
    riva_material_parameters_t material = {};
    material.h=p->h; material.m=p->m; material.M=p->M;
    material.kd=p->kd; material.zeta=p->zeta;
    material.e_max=p->e_max; material.e_min=p->e_min;
    material.Q=p->Q; material.R=p->R; material.n_G=p->n_G;
    return material;
}

RIVA_HD static inline int riva_material_parameters_valid(
    const riva_parameters_t *p,const riva_material_parameters_t *material)
{
    return material && isfinite(material->h) && material->h>0.0 &&
        isfinite(material->m) && material->m>=0.0 &&
        isfinite(material->M) && material->M>0.0 &&
        isfinite(material->kd) && material->kd>0.0 &&
        isfinite(material->zeta) && material->zeta>=0.0 &&
        isfinite(material->e_max) && isfinite(material->e_min) &&
        material->e_max>material->e_min && material->e_min>=0.0 &&
        p->state_shakedown_reference_void_ratio>=material->e_min &&
        p->state_shakedown_reference_void_ratio<=material->e_max &&
        isfinite(material->Q) && material->Q>0.0 &&
        isfinite(material->R) && material->R>=0.0 &&
        isfinite(material->n_G) && material->n_G>=0.0 &&
        isfinite(p->p_residual) && p->p_residual>=0.0 &&
        (p->geostatic_admission_enabled==0 ||
         p->geostatic_admission_enabled==1) &&
        (p->v11_enabled==0 ||
         (isfinite(p->v11_k) && p->v11_k>=0.0 &&
          isfinite(p->v11_eta_floor) && p->v11_eta_floor>0.0 &&
          p->v11_eta_floor<=1.0)) &&
        material->n_G<=1.0;
}

/* Relative density is the external state descriptor for both the frozen and
 * custom V8 interfaces. Keep void ratio as an internal constitutive variable
 * and perform this affine conversion at the input boundary. */
RIVA_HD static inline double riva_void_ratio_from_material_relative_density(
    const riva_material_parameters_t *material,double relative_density)
{
    return material->e_max-relative_density*(material->e_max-material->e_min);
}

RIVA_HD static inline double riva_void_ratio_from_relative_density(
    const riva_parameters_t *p, double relative_density)
{
    const riva_material_parameters_t material=
        riva_reference_material_parameters(p);
    return riva_void_ratio_from_material_relative_density(
        &material,relative_density);
}

RIVA_HD static inline double riva_relative_state(
    const riva_parameters_t *p,const riva_material_parameters_t *material,
    double pressure,double void_ratio)
{
    if (!p->state_enabled) return 0.0;
    const double dr=riva_clip((material->e_max-void_ratio)/
        (material->e_max-material->e_min),0.0,1.0);
    const double arg=riva_max(100.0*riva_max(pressure,riva_cone_pressure_floor(p))/p->p_atm,1.0e-12);
    const double critical=riva_clip(material->R/
        riva_max(material->Q-log(arg),0.10),0.0,1.0);
    return riva_clip(dr-critical,-p->relative_state_limit,p->relative_state_limit);
}

RIVA_HD static inline void riva_surfaces(const riva_parameters_t *p,
    const riva_material_parameters_t *material,
    double pressure,double void_ratio,double *mb,double *md,double *xi)
{
    *xi=riva_relative_state(p,material,pressure,void_ratio);
    *mb=material->M*exp(p->n_b*(*xi));
    *md=material->kd*exp(-p->n_d*(*xi));
}

/* Keep a geostatically admitted over-bound stress non-expansive while the
 * translated cone shrinks: it may remain outside the calibrated cone, but its
 * normalized distance from the apex may not grow.  This is a radial, closed-
 * form correction and therefore adds no local iterations or substeps. */
RIVA_HD static inline riva_tensor_t riva_regularize_deviator(
    const riva_parameters_t *p,const riva_material_parameters_t *material,
    double pressure,double void_ratio,riva_tensor_t prior_alpha,
    riva_tensor_t deviator,riva_tensor_t *alpha_out)
{
    double mb,md,xi;
    riva_surfaces(p,material,pressure,void_ratio,&mb,&md,&xi);
    (void)md; (void)xi;
    riva_tensor_t alpha=riva_scale(deviator,1.0/pressure);
    double limit=sqrt(2.0/3.0)*mb;
    const int admit_overbound=p->geostatic_admission_enabled||
        p->admit_inherited_overbound;
    if (admit_overbound)
        limit=riva_max(limit,riva_norm(prior_alpha));
    const double norm=riva_norm(alpha);
    if (norm>limit) {
        double allowed=limit;
        if (admit_overbound) {
            /* C1-smooth cap: value- and slope-continuous at the limit,
             * overshoot bounded by the band width. The hard snap leaves a
             * derivative kink exactly where admitted over-bound states sit
             * after their first commit (the raised limit equals the
             * committed ratio), and implicit equilibrium chatters on that
             * kink. Gated: the frozen path keeps the hard snap bit-exact. */
            const double band=0.02*limit;
            const double excess=norm-limit;
            allowed=limit+band*(1.0-exp(-excess/band));
        }
        alpha=riva_scale(alpha,allowed/norm);
        deviator=riva_scale(alpha,pressure);
    }
    if (alpha_out) *alpha_out=alpha;
    return deviator;
}

RIVA_HD static inline void riva_moduli(const riva_parameters_t *p,
    const riva_material_parameters_t *material,double pressure,
    double *shear,double *bulk)
{
    const double ratio=riva_max(pressure,riva_cone_pressure_floor(p))/p->p_ref;
    const double young=p->E_ref*pow(ratio,material->n_G);
    *shear=young/(2.0*(1.0+p->nu));
    *bulk=young/(3.0*(1.0-2.0*p->nu));
}

RIVA_HD static inline double riva_state_factor(const riva_parameters_t *p,double xi)
{
    if (!p->state_contraction_enabled) return 1.0;
    const double d=xi-p->reference_relative_state;
    const double exponent=d>=0.0?p->dense_state_contraction_exponent:p->loose_state_contraction_exponent;
    return riva_clip(exp(-exponent*d),p->state_factor_minimum,p->state_factor_maximum);
}

RIVA_HD static inline double riva_effective_knee(const riva_parameters_t *p,double xi)
{
    if (!p->state_dependent_knee_enabled) return p->low_amplitude_knee_ratio;
    const double d=xi-p->reference_relative_state;
    const double exponent=d>=0.0?p->dense_knee_state_exponent:p->loose_knee_state_exponent;
    return riva_clip(p->low_amplitude_knee_ratio*exp(exponent*d),p->effective_knee_minimum,1.0);
}

RIVA_HD static inline double riva_amplitude_factor(const riva_parameters_t *p,
    double amplitude,double knee)
{
    if (!p->amplitude_scaling_enabled) return 1.0;
    const double ratio=riva_max(amplitude/p->cyclic_amplitude_reference,1.0e-12);
    const double f=pow(ratio,p->cyclic_amplitude_exponent)*
        pow(riva_min(1.0,ratio/knee),p->low_amplitude_exponent);
    return riva_clip(f,p->amplitude_factor_minimum,p->amplitude_factor_maximum);
}

RIVA_HD static inline double riva_projected_bias(const riva_state_t *s)
{
    if (!s->cyclic_phase_active) return 0.0;
    const double norm=riva_norm(s->cyclic_direction);
    return norm<=1.0e-12?0.0:fabs(riva_ddot(s->static_bias_tensor,s->cyclic_direction))/norm;
}

RIVA_HD static inline double riva_initial_relative_density(
    const riva_parameters_t *p,const riva_material_parameters_t *material,
    const riva_state_t *s)
{
    const double arg=riva_max(100.0*riva_max(s->pressure_anchor,riva_cone_pressure_floor(p))/p->p_atm,1.0e-12);
    return riva_clip(s->initial_relative_state+material->R/
        riva_max(material->Q-log(arg),0.10),0.0,1.0);
}

RIVA_HD static inline double riva_reference_relative_density(
    const riva_parameters_t *p,const riva_material_parameters_t *material)
{ return (material->e_max-p->state_shakedown_reference_void_ratio)/
    (material->e_max-material->e_min); }

RIVA_HD static inline double riva_raw_shakedown(const riva_parameters_t *p,
    const riva_material_parameters_t *material,double pressure,
    const riva_state_t *s)
{
    if (!p->state_shakedown_enabled || !s->cyclic_phase_active ||
        s->amplitude_reversals<p->state_shakedown_minimum_reversals) return 1.0;
    const double bias=riva_projected_bias(s);
    if (bias<=1.0e-14) return 1.0;
    const double distance=riva_initial_relative_density(p,material,s)-
        riva_reference_relative_density(p,material);
    const double state_term=exp(p->state_shakedown_state_sensitivity*distance);
    const double anchor_release=1.0-exp(-pow(distance/p->state_shakedown_anchor_width,2.0));
    const double bias_term=pow(bias/p->state_shakedown_bias_reference,
                               p->state_shakedown_bias_exponent);
    const double pressure_term=pow(riva_clip(pressure/riva_max(s->pressure_anchor,riva_cone_pressure_floor(p)),0.0,1.0),
                                   p->state_shakedown_pressure_exponent);
    return riva_clip(1.0+p->state_shakedown_scale*anchor_release*state_term*
                    bias_term*pressure_term,1.0,p->state_shakedown_factor_maximum);
}

RIVA_HD static inline double riva_dense_state_weight(const riva_parameters_t *p,
    const riva_material_parameters_t *material,const riva_state_t *s)
{
    const double d=riva_max(riva_initial_relative_density(p,material,s)-
        riva_reference_relative_density(p,material),0.0);
    return 1.0-exp(-pow(d/p->state_shakedown_anchor_width,2.0));
}

RIVA_HD static inline double riva_shakedown(const riva_parameters_t *p,
    const riva_material_parameters_t *material,double pressure,
    const riva_state_t *s)
{
    const double raw=riva_raw_shakedown(p,material,pressure,s);
    const double excess=riva_max(s->cyclic_amplitude-
        p->state_shakedown_dense_hardening_amplitude_onset,0.0);
    const double dense_scale=p->state_shakedown_dense_hardening_scale*
        exp(-p->state_shakedown_dense_hardening_amplitude_decay*excess);
    const double branch=1.0+riva_dense_state_weight(p,material,s)*(dense_scale-1.0);
    return riva_clip(1.0+(raw-1.0)*branch,1.0,p->state_shakedown_factor_maximum);
}

RIVA_HD static inline double riva_shakedown_compliance(const riva_parameters_t *p,
    const riva_material_parameters_t *material,double pressure,
    const riva_state_t *s)
{
    const double raw=riva_raw_shakedown(p,material,pressure,s);
    const double d=riva_max(riva_initial_relative_density(p,material,s)-
        riva_reference_relative_density(p,material),0.0);
    const double attenuation=exp(-p->state_shakedown_dense_compliance_decay*d*d);
    return 1.0+(raw-1.0)*attenuation;
}

RIVA_HD static inline void riva_cyclic_flow(const riva_parameters_t *p,
    double pressure,const riva_state_t *s,double *shear_factor,double *hardening_factor)
{
    if (!p->cyclic_flow_correction_enabled || p->diag_no_cyclic_flow ||
        p->v11_enabled ||
        s->amplitude_reversals<p->cyclic_flow_minimum_reversals) {
        *shear_factor=1.0; *hardening_factor=1.0; return;
    }
    /* NOTE (V9.2 attempts, 2026-08-25): scaling this cut by the half-cycle
       plastic activity fixes the moderate-amplitude damping hole (xi 0.5%
       -> 11.7% at gamma = 0.1%) and the small-strain over-softness
       (G_sec/G_el 0.14 -> 0.4), BUT regresses either the stress-controlled
       cyclic-strength curve (ep_ref = 1.2e-3: CSR 0.30 stops triggering) or
       the biased ratcheting direction (ep_ref = 3e-4: 3-9% UPSLOPE creep at
       alpha >= 0.25). The cut is load-bearing for the calibrated response;
       recalibrating it needs lab anchors (G/gamma + damping + K_alpha
       data) -- model-owner territory. Kept frozen. */
    const double ratio=riva_clip(pressure/riva_max(s->pressure_anchor,riva_cone_pressure_floor(p)),0.0,1.0);
    const double activity=pow(ratio,p->cyclic_flow_pressure_exponent);
    double w_cut=1.0;
    if (p->l3_ep_ref>0.0) {
        const double act=riva_min(1.0,s->ep_half_last/p->l3_ep_ref);
        w_cut=p->l3_cut_floor+(1.0-p->l3_cut_floor)*act;
    }
    *shear_factor=1.0-p->cyclic_shear_modulus_reduction*activity*w_cut;
    *hardening_factor=1.0+p->cyclic_hardening_boost*activity;
}

RIVA_HD static inline double riva_bias_hardening(const riva_parameters_t *p,
    const riva_state_t *s)
{
    if (!p->static_bias_enabled || !s->cyclic_phase_active ||
        s->amplitude_reversals<p->bias_minimum_reversals) return 0.0;
    const double bias=riva_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    const double amplitude=p->bias_amplitude_ratio*s->cyclic_amplitude;
    const double margin=riva_max(bias-amplitude,0.0);
    const double crossing=riva_max(amplitude/bias-1.0,0.0);
    return p->bias_hardening_intercept*pow(bias,p->bias_intercept_exponent)*
        exp(-p->bias_crossing_decay*crossing)+
        p->bias_hardening_scale*pow(margin,p->bias_margin_exponent);
}

RIVA_HD static inline void riva_moduli_for_state(const riva_parameters_t *p,
    const riva_material_parameters_t *material,double pressure,
    const riva_state_t *s,double *shear,double *bulk)
{
    double sf,hf; riva_moduli(p,material,pressure,shear,bulk);
    riva_cyclic_flow(p,pressure,s,&sf,&hf); (void)hf;
    *shear *= sf;
    const double compliance=riva_shakedown_compliance(p,material,pressure,s);
    *shear *= riva_clip(pow(compliance,-p->state_shakedown_compliance_exponent),
                       p->state_shakedown_shear_multiplier_minimum,1.0);
}

RIVA_HD static inline double riva_hardening_for_state(const riva_parameters_t *p,
    const riva_material_parameters_t *material,double pressure,
    const riva_state_t *s)
{
    const double ratio=riva_max(pressure,riva_cone_pressure_floor(p))/p->p_ref;
    double value=material->h*pow(ratio,-p->q_H);
    double sf,hf; riva_cyclic_flow(p,pressure,s,&sf,&hf); (void)sf;
    value *= hf;
    double boost=p->diag_no_bias_hardening?0.0:riva_bias_hardening(p,s);
    if (p->l3_ep_ref>0.0) {
        const double act=riva_min(1.0,s->ep_half_last/p->l3_ep_ref);
        boost *= p->l3_boost_floor+(1.0-p->l3_boost_floor)*act;
    }
    if (p->v11_enabled) {
        double mb_m,md_m,xi_m;
        riva_surfaces(p,material,pressure,s->void_ratio,&mb_m,&md_m,&xi_m);
        (void)md_m; (void)xi_m;
        const double eta_ratio=riva_clip(
            s->eta_m/riva_max(mb_m,1.0e-6),p->v11_eta_floor,1.0);
        value *= pow(1.0/eta_ratio,p->v11_k);
        boost=0.0;   /* memory factor replaces the bias boost */
    }
    const double pressure_ratio=riva_clip(pressure/riva_max(s->pressure_anchor,riva_cone_pressure_floor(p)),0.0,1.0);
    const double confinement=pow(riva_clip(p->bias_reference_pressure/
        riva_max(s->pressure_anchor,riva_cone_pressure_floor(p)),0.25,4.0),p->bias_confinement_exponent);
    value *= 1.0+boost*confinement*pow(pressure_ratio,p->bias_pressure_exponent);
    return value*riva_shakedown(p,material,pressure,s);
}

RIVA_HD static inline double riva_bulk_reference(const riva_parameters_t *p,
    const riva_material_parameters_t *material)
{ double g,k; riva_moduli(p,material,p->p_ref,&g,&k); return k; }

RIVA_HD static inline double riva_confining_strain(const riva_parameters_t *p,
    const riva_material_parameters_t *material,double pressure,double anchor)
{
    const double pressure_floor=riva_cone_pressure_floor(p);
    pressure=riva_max(pressure,pressure_floor); anchor=riva_max(anchor,pressure_floor);
    const double bulk=riva_bulk_reference(p,material);
    if (fabs(material->n_G-1.0)<1.0e-12)
        return p->p_ref/bulk*log(anchor/pressure);
    const double factor=pow(p->p_ref,material->n_G)/
        (bulk*(1.0-material->n_G));
    return factor*(pow(anchor,1.0-material->n_G)-
        pow(pressure,1.0-material->n_G));
}

RIVA_HD static inline double riva_pressure_from_confining(
    const riva_parameters_t *p,const riva_material_parameters_t *material,
    double eps_vc,double anchor,int *at_floor)
{
    const double pressure_floor=riva_cone_pressure_floor(p);
    const double threshold=riva_confining_strain(p,material,pressure_floor,anchor);
    if (eps_vc>=threshold) { *at_floor=1; return pressure_floor; }
    *at_floor=0;
    const double bulk=riva_bulk_reference(p,material);
    double pressure;
    if (fabs(material->n_G-1.0)<1.0e-12)
        pressure=anchor*exp(-bulk*eps_vc/p->p_ref);
    else {
        double base=pow(anchor,1.0-material->n_G)-
            (1.0-material->n_G)*bulk*eps_vc/
                pow(p->p_ref,material->n_G);
        base=riva_max(base,pow(pressure_floor,1.0-material->n_G));
        pressure=pow(base,1.0/(1.0-material->n_G));
    }
    return riva_max(pressure,pressure_floor);
}

RIVA_HD static inline int riva_initialize_material(const riva_parameters_t *p,
    const riva_material_parameters_t *material,riva_tensor_t stress,
    double void_ratio,riva_state_t *state)
{
    const double pressure=riva_cone_pressure(p,stress);
    if (!riva_material_parameters_valid(p,material) || !(pressure>0.0) ||
        !isfinite(pressure) || !isfinite(void_ratio)) return 0;
    *state=riva_state_t{};
    const riva_tensor_t dev=riva_dev(stress);
    double mb,md,xi;
    riva_surfaces(p,material,pressure,void_ratio,&mb,&md,&xi); (void)mb;
    state->stress=stress; state->alpha=riva_scale(dev,1.0/pressure);
    state->D_ir=p->irreversible_enabled?
        material->zeta*sqrt(2.0/3.0)*md:0.0;
    state->D=state->D_ir; state->beta=p->beta0;
    state->void_ratio=void_ratio; state->pressure_anchor=pressure;
    state->last_reversal_deviator=dev; state->amplitude_factor=1.0;
    state->initial_relative_state=xi; state->state_contraction_factor=riva_state_factor(p,xi);
    state->effective_knee_ratio=riva_effective_knee(p,xi);
    state->geostatic_deviator=dev; state->initialized=1;
    state->eta_m=sqrt(1.5)*riva_norm(state->alpha);
    return 1;
}

RIVA_HD static inline int riva_initialize(const riva_parameters_t *p,
    riva_tensor_t stress,double void_ratio,riva_state_t *state)
{
    const riva_material_parameters_t material=
        riva_reference_material_parameters(p);
    return riva_initialize_material(p,&material,stress,void_ratio,state);
}

RIVA_HD static inline int riva_begin_dynamic_phase(const riva_parameters_t *p,
    riva_tensor_t reference_stress,riva_state_t *state)
{
    if (!state->initialized) return 0;
    const riva_tensor_t current=riva_dev(state->stress);
    const riva_tensor_t reference=riva_dev(reference_stress);
    const riva_tensor_t bias=riva_sub(current,reference);
    const double anchor=riva_max(state->pressure_anchor,riva_cone_pressure_floor(p));
    state->static_bias_tensor=riva_scale(bias,1.0/anchor);
    state->static_bias_index=riva_norm(bias)/anchor;
    state->cyclic_direction=riva_zero(); state->cyclic_phase_active=1;
    state->bias_ratchet_strain=0.0; state->bias_reversible_volume=0.0;
    state->eps_v_total=state->physical_eps_v_total;
    state->last_host_deviatoric_strain_direction=riva_zero();
    return 1;
}

RIVA_HD static inline void riva_dilatancy(const riva_parameters_t *p,
    const riva_material_parameters_t *material,
    riva_tensor_t alpha,riva_tensor_t normal,double beta,riva_tensor_t fabric,
    double pressure,double void_ratio,double eps_v_ir,double eps_v_re,
    double *D_ir,double *D_re)
{
    double mb,md,xi;
    riva_surfaces(p,material,pressure,void_ratio,&mb,&md,&xi); (void)xi;
    const double eta=riva_ddot(alpha,normal);
    const double eta_d=sqrt(2.0/3.0)*md;
    *D_ir=0.0;
    if (p->irreversible_enabled) {
        const double distance=riva_max(eta_d-eta,0.0);
        const double amplitude=material->zeta*
            (1.0+riva_max(riva_ddot(fabric,normal),0.0));
        double gate=1.0;
        if (p->contraction_gate_enabled && p->C_D>0.0) {
            const double reversal=(sqrt(2.0/3.0)*mb-eta)/riva_max(beta,1.0e-12);
            const double c_in=p->C_in_ratio*sqrt(2.0/3.0)*mb;
            gate=(reversal+c_in)*(reversal+c_in)/(reversal*reversal+p->C_D);
        }
        *D_ir=amplitude*gate*distance*
            exp(-p->irreversible_decay*riva_max(-eps_v_ir,0.0));
    }
    *D_re=0.0;
    if (p->reversible_enabled) {
        if (eta>eta_d)
            *D_re=-material->zeta*p->reversible_ratio*(eta-eta_d);
        else if (eps_v_re>0.0)
            *D_re=riva_min(p->reversible_release_rate*eps_v_re,riva_max(eta_d,1.0e-12));
    }
}

RIVA_HD static inline double riva_irreversible_factor(const riva_parameters_t *p,
    const riva_state_t *s)
{
    if (p->v12_enabled) {
        /* pure state-coupled chain - no amplitude factor, no bias boost */
        double factor=s->state_contraction_factor*p->v12_dir_scale;
        if (p->v12_chi>0.0 && s->eps_v_irreversible<-1.0e-12) {
            const double chi=riva_min(1.0,
                p->v12_chi*riva_max(s->eps_v_reversible,0.0)/
                fabs(s->eps_v_irreversible));
            factor *= exp(chi);
        }
        if (p->v11_rdr>0.0) {
            const double w=p->v11_rdr/
                (p->v11_rdr+riva_max(s->ep_eq_since_reversal,0.0));
            factor *= w*w;
        }
        if (p->v11_eta_ir>0.0)
            factor *= exp(-p->v11_eta_ir*fabs(s->eps_v_irreversible));
        return factor;
    }
    double amp=s->amplitude_factor;
    if (p->v11_enabled && amp<p->v11_dir_floor) amp=p->v11_dir_floor;
    double factor=amp*s->state_contraction_factor;
    const double bias=riva_projected_bias(s);
    factor *= 1.0+p->bias_contraction_scale*pow(bias,p->bias_contraction_exponent);
    if (p->v11_enabled) {
        if (p->v11_rdr>0.0) {
            const double w=p->v11_rdr/
                (p->v11_rdr+riva_max(s->ep_eq_since_reversal,0.0));
            factor *= w*w;
        }
        if (p->v11_eta_ir>0.0)
            factor *= exp(-p->v11_eta_ir*fabs(s->eps_v_irreversible));
    }
    return factor;
}

RIVA_HD static inline void riva_register_reversal(const riva_parameters_t *p,
    const riva_state_t *old,riva_tensor_t *alpha0,riva_tensor_t *alpha01,
    double *ep_eq,double *beta,int64_t *reversals)
{
    const double ratio=(*ep_eq)/p->overshoot_strain;
    const double modifier=ratio>1.0?0.0:1.0-ratio;
    const riva_tensor_t previous=*alpha0;
    *alpha0=riva_norm(previous)==0.0?old->alpha:
        riva_add(riva_scale(*alpha01,modifier),riva_scale(old->alpha,1.0-modifier));
    *alpha01=previous; *ep_eq=0.0; *beta=p->beta0;
    *reversals=old->reversals+1;
}

RIVA_HD static inline riva_state_t riva_backbone_forward_euler(
    const riva_parameters_t *p,const riva_material_parameters_t *material,
    const riva_state_t *old,riva_tensor_t deps,int force_reversal,
    int allow_legacy_reversal)
{
    riva_state_t state=*old;
    const riva_tensor_t s_t=riva_dev(old->stress);
    const double pressure_floor=riva_cone_pressure_floor(p);
    const double p_t=riva_max(riva_cone_pressure(p,old->stress),pressure_floor);
    double shear,bulk;
    riva_moduli_for_state(p,material,p_t,old,&shear,&bulk);
    const double h_eff=riva_hardening_for_state(p,material,p_t,old);
    const double threshold=riva_confining_strain(
        p,material,pressure_floor,old->pressure_anchor);
    const int floor_active=p->compatibility_enabled && old->eps_v_confining>=threshold;
    const double coupling_bulk=floor_active?0.0:bulk;

    const double deps_v=riva_trace(deps);
    const riva_tensor_t deps_d=riva_sub(deps,riva_iso(deps_v/3.0));
    const riva_tensor_t s_trial=riva_add(s_t,riva_scale(deps_d,2.0*shear));
    double p_trial;
    if (p->compatibility_enabled) {
        int at_floor=0;
        p_trial=riva_pressure_from_confining(
            p,material,old->eps_v_confining+deps_v,
            old->pressure_anchor,&at_floor);
    } else p_trial=riva_max(p_t-bulk*deps_v,pressure_floor);
    if (p->p_max>0.0 && p_trial>p->p_max) p_trial=p->p_max;
    const riva_tensor_t alpha_trial=riva_scale(s_trial,1.0/riva_max(p_trial,pressure_floor));
    const riva_tensor_t delta_alpha_trial=riva_sub(alpha_trial,old->alpha);

    riva_tensor_t alpha0=old->alpha0,alpha01=old->alpha01;
    double beta_t=old->beta,ep_eq=old->ep_eq_since_reversal;
    int64_t reversals=old->reversals;
    const double legacy=riva_ddot(riva_sub(alpha_trial,alpha0),delta_alpha_trial);
    const int reversal=force_reversal || (allow_legacy_reversal && legacy<0.0);
    if (reversal) riva_register_reversal(p,old,&alpha0,&alpha01,&ep_eq,&beta_t,&reversals);

    const double alpha_n_t=riva_ddot(old->alpha,old->n);
    double beta_h=riva_max(beta_t,riva_max(p->beta_floor,1.0e-12));
    if (p->beta_cap>0.0 && beta_h>p->beta_cap) beta_h=p->beta_cap;
    double hardening_beta=pow(beta_h,material->m);
    if (p->beta_reserve>0.0 && old->reversals==0) {
        /* State-dependent gate: the reserve engages only when the committed
         * stress ratio sits within 80% of the bounding surface ??the regime
         * where beta (and with it the plastic modulus) collapses. K0
         * consolidation states (~70%) never engage it; wave-driven points
         * engage it exactly when driven near the bound. Cyclic stages run
         * with the reserve parameter zeroed by the host regardless. */
        double mb_r,md_r,xi_r;
        riva_surfaces(p,material,p_t,old->void_ratio,&mb_r,&md_r,&xi_r);
        if (riva_norm(old->alpha)>0.8*sqrt(2.0/3.0)*mb_r)
            hardening_beta+=pow(p->beta_reserve,material->m);
    }
    const double hardening=p_t*h_eff*hardening_beta;
    double denominator=2.0*shear+(2.0/3.0)*hardening-
        coupling_bulk*old->D*alpha_n_t;
    const double floor=p->denominator_floor_ratio*2.0*shear;
    int64_t denominator_hits=old->denominator_floor_hits;
    if (denominator<floor) { denominator=floor; denominator_hits++; }
    const double numerator=2.0*shear*riva_ddot(deps_d,old->n)+
        coupling_bulk*deps_v*alpha_n_t;
    const double d_lambda=riva_max(0.0,numerator/denominator);
    riva_tensor_t s=riva_add(s_t,riva_scale(riva_sub(deps_d,riva_scale(old->n,d_lambda)),
                                     2.0*shear));

    const double d_ir=riva_max(old->D_ir,0.0);
    double d_re=old->D_re;
    if (d_lambda>0.0 && d_re>0.0)
        d_re=riva_min(d_re,old->eps_v_reversible/d_lambda);
    const double eps_v_ir=old->eps_v_irreversible-d_lambda*d_ir;
    const double eps_v_re=riva_max(old->eps_v_reversible-d_lambda*d_re,0.0);
    const double eps_v_total=old->eps_v_total+deps_v;
    const double eps_v_confining=eps_v_total-eps_v_ir-eps_v_re;
    double pressure; int at_floor=0;
    if (p->compatibility_enabled)
        pressure=riva_pressure_from_confining(
            p,material,eps_v_confining,old->pressure_anchor,&at_floor);
    else {
        pressure=riva_max(p_t-bulk*(deps_v+d_lambda*(d_ir+d_re)),pressure_floor);
        at_floor=pressure<=pressure_floor;
    }
    if (p->p_max>0.0 && pressure>p->p_max) pressure=p->p_max;
    riva_tensor_t stress=riva_sub(
        s,riva_iso(riva_physical_pressure(p,pressure)));
    ep_eq += d_lambda;
    const double void_ratio=old->void_ratio+(1.0+old->void_ratio)*deps_v;
    riva_tensor_t fabric=old->fabric;
    if (p->fabric_enabled && d_lambda>0.0) {
        const riva_tensor_t term=riva_add(riva_scale(old->n,sqrt(2.0/3.0)*p->z_max),fabric);
        fabric=riva_add(fabric,riva_scale(term,-p->c_z*d_lambda*riva_max(-d_re,0.0)));
    }

    double mb,md,xi;
    riva_surfaces(p,material,pressure,void_ratio,&mb,&md,&xi);
    (void)md; (void)xi;
    riva_tensor_t alpha;
    s=riva_regularize_deviator(
        p,material,pressure,void_ratio,old->alpha,s,&alpha);
    stress=riva_sub(s,riva_iso(riva_physical_pressure(p,pressure)));
    const riva_tensor_t offset=riva_sub(alpha,alpha0);
    const double qa=riva_ddot(offset,offset);
    const double qb=2.0*riva_ddot(offset,alpha);
    const double qc=riva_ddot(alpha,alpha)-(2.0/3.0)*mb*mb;
    const double disc=qb*qb-4.0*qa*qc;
    int64_t beta_fallbacks=old->beta_fallbacks;
    double beta;
    if (qa<=1.0e-14 || disc<0.0) { beta=p->beta0; beta_fallbacks++; }
    else {
        beta=(-qb+sqrt(riva_max(disc,0.0)))/(2.0*qa);
        if (!isfinite(beta)) { beta=p->beta0; beta_fallbacks++; }
    }
    beta=riva_max(beta,1.0e-6);
    const riva_tensor_t mapped=riva_add(alpha,riva_scale(offset,beta));
    const double mapped_norm=riva_norm(mapped);
    const riva_tensor_t normal=mapped_norm>1.0e-12?riva_scale(mapped,1.0/mapped_norm):riva_zero();
    double D_ir,D_re;
    riva_dilatancy(p,material,alpha,normal,beta,fabric,pressure,void_ratio,
                  eps_v_ir,eps_v_re,&D_ir,&D_re);
    D_ir *= riva_irreversible_factor(p,old);

    state.stress=stress; state.alpha=alpha; state.alpha0=alpha0;
    state.alpha01=alpha01; state.n=normal; state.fabric=fabric;
    state.eta_m=riva_max(old->eta_m,sqrt(1.5)*riva_norm(alpha));
    state.D_ir=D_ir; state.D_re=D_re; state.D=D_ir+D_re; state.beta=beta;
    state.lambda_total=old->lambda_total+d_lambda;
    state.ep_eq_since_reversal=ep_eq; state.void_ratio=void_ratio;
    state.reversals=reversals;
    state.pressure_floor_hits=old->pressure_floor_hits+(at_floor?1:0);
    state.denominator_floor_hits=denominator_hits;
    state.beta_fallbacks=beta_fallbacks;
    state.eps_v_total=eps_v_total; state.eps_v_confining=eps_v_confining;
    state.eps_v_irreversible=eps_v_ir; state.eps_v_reversible=eps_v_re;
    if (state.reversals>old->reversals) {
        const riva_tensor_t reversal_dev=riva_dev(old->stress);
        const double excursion=riva_norm(riva_sub(reversal_dev,old->last_reversal_deviator));
        const double divisor=old->amplitude_reversals==0?1.0:2.0;
        const double amplitude=excursion/(divisor*riva_max(old->pressure_anchor,pressure_floor));
        if (amplitude>1.0e-12) {
            state.cyclic_amplitude=amplitude;
            state.amplitude_factor=riva_amplitude_factor(p,amplitude,old->effective_knee_ratio);
        }
        state.last_reversal_deviator=reversal_dev;
        state.amplitude_reversals=old->amplitude_reversals+1;
        /* plastic activity of the half-cycle that just completed: the
           amplitude-consistency discriminant for the bias reversible-volume
           block (PM4Sand/SANISAND-style: cyclic volumetric response follows
           the plastic multiplier, so micro/elastic cycles drive nothing) */
        state.ep_half_last=old->ep_eq_since_reversal;
    }
    return state;
}

RIVA_HD static inline double riva_bias_reversible_factor(
    const riva_parameters_t *p,const riva_state_t *s)
{
    if (!p->static_bias_enabled || !s->cyclic_phase_active) return 1.0;
    return 1.0+p->bias_reversible_scale*
        pow(riva_projected_bias(s),p->bias_reversible_exponent);
}

RIVA_HD static inline double riva_smoothstep(double value)
{ value=riva_clip(value,0.0,1.0); return value*value*(3.0-2.0*value); }

RIVA_HD static inline double riva_bias_reversible_volume_target(
    const riva_parameters_t *p,const riva_state_t *s)
{
    if (!p->static_bias_enabled || !p->bias_reversible_volume_enabled ||
        !s->cyclic_phase_active || s->amplitude_reversals<1) return 0.0;
    const double direction_norm=riva_norm(s->cyclic_direction);
    if (direction_norm<=1.0e-14 || s->cyclic_amplitude<=1.0e-14) return 0.0;
    const double bias=riva_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    /* Plastic-activity gate (V9): the response follows the plastic
       multiplier of the last completed half-cycle, PM4Sand/SANISAND-style ??
       elastic/micro cycles drive nothing, structurally, with no dependence
       on the current pressure (no p-collapse feedback loop). ep_ref sits
       below the golden floor (3.4e-5 in cyclic_bias0375_dense), so every
       calibrated evaluation saturates at exactly 1.0. */
    /* V9: response proportional to the last half-cycle's plastic activity */
    const double ep_scale=riva_min(1.0,
        s->ep_half_last/p->bias_reversible_volume_ep_cal);
    if (ep_scale<=0.0) return 0.0;
    const riva_tensor_t phase_direction=riva_scale(s->cyclic_direction,1.0/direction_norm);
    const riva_tensor_t current=riva_dev(s->stress);
    const riva_tensor_t static_dev=riva_add(s->geostatic_deviator,
        riva_scale(s->static_bias_tensor,s->pressure_anchor));
    const double dynamic_projection=riva_ddot(riva_sub(current,static_dev),phase_direction);
    const double phase=riva_clip(dynamic_projection/
        (s->pressure_anchor*s->cyclic_amplitude),-1.0,1.0);
    /* V9 knee: frozen expression at/above the knee, linear fade below.
     * NOTE (integration, 2026-08-25): upstream research/rivasand-geostatic
     * instead fades the raw inverse-bias law with a smoothstep^3 below
     * 0.25*reference_bias. The V9 knee subsumes it: it also bounds the
     * factor (cap ~4.6 at the 0.30 knee) through the 0.16-0.30 band where
     * the raw law still reaches ~16x, and it is what the re-frozen golden
     * oracle and the P1 element validation were run against. Flagged for
     * the PR discussion rather than silently double-gating. */
    const double bias_eff=riva_max(bias,p->bias_reversible_volume_bias_knee);
    const double bias_factor=pow(p->bias_reversible_volume_reference_bias/bias_eff,
                                 p->bias_reversible_volume_bias_exponent)*
        riva_min(1.0,bias/p->bias_reversible_volume_bias_knee);
    const double pressure_factor=pow(riva_max(s->pressure_anchor,riva_cone_pressure_floor(p))/
        p->bias_reference_pressure,p->bias_reversible_volume_pressure_exponent);
    const double oscillation_buildup=1.0-exp(-(double)s->amplitude_reversals/
        p->bias_reversible_volume_buildup_reversals);
    const double oscillation=-p->bias_reversible_volume_amplitude*bias_factor*
        pressure_factor*oscillation_buildup*phase;
    const double mean_buildup=1.0-exp(-(double)s->amplitude_reversals/
        p->bias_reversible_mean_buildup_reversals);
    const double mean_shift=p->bias_reversible_mean_scale*
        (p->bias_reversible_mean_transition_pressure-s->pressure_anchor)/
        p->bias_reference_pressure*(bias/p->bias_reversible_volume_reference_bias)*
        mean_buildup;
    /* low-pressure fade: exactly 1 at the calibrated anchors (>= fade
       pressure), tapering the fixed volumetric amplitude where the small
       bulk modulus would otherwise turn it into confinement-scale swings */
    const double lowp_fade=riva_smoothstep(
        riva_max(s->pressure_anchor,riva_cone_pressure_floor(p))/
        p->bias_reversible_volume_pressure_fade);
    return (oscillation+mean_shift)*ep_scale*lowp_fade;
}

RIVA_HD static inline double riva_bias_ratchet_activity(
    const riva_parameters_t *p,const riva_state_t *s)
{
    if (!p->static_bias_enabled || !p->bias_ratchet_enabled ||
        !s->cyclic_phase_active || s->amplitude_reversals<p->bias_minimum_reversals)
        return 0.0;
    const double bias=riva_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    const double amplitude_position=(s->cyclic_amplitude-p->bias_ratchet_amplitude_onset)/
        (p->bias_ratchet_amplitude_full-p->bias_ratchet_amplitude_onset);
    const double ratio_position=(s->cyclic_amplitude/bias-p->bias_ratchet_ratio_full)/
        (p->bias_ratchet_ratio_cutoff-p->bias_ratchet_ratio_full);
    return riva_smoothstep(amplitude_position)*(1.0-riva_smoothstep(ratio_position));
}

RIVA_HD static inline double riva_bias_ratchet_capacity(
    const riva_parameters_t *p,const riva_state_t *s)
{
    const double bias=riva_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    const double factor=pow(p->bias_ratchet_reference_bias/bias,
                            p->bias_ratchet_bias_exponent);
    /* The repeated factor is intentional and preserves the frozen V8 oracle. */
    return p->bias_ratchet_limit*riva_bias_ratchet_activity(p,s)*factor*factor;
}

RIVA_HD static inline double riva_bias_ratchet_rate_factor(
    const riva_parameters_t *p,const riva_state_t *s)
{
    const double capacity=riva_bias_ratchet_capacity(p,s);
    if (capacity<=1.0e-14 || s->bias_ratchet_strain>=capacity) return 0.0;
    const double stress_level=pow(riva_max(s->pressure_anchor/
        p->bias_reference_pressure-1.0,0.0),p->bias_ratchet_pressure_exponent);
    const double bias=riva_projected_bias(s);
    const double bias_factor=pow(p->bias_ratchet_reference_bias/
        riva_max(bias,1.0e-14),p->bias_ratchet_bias_exponent);
    const double saturation=riva_max(1.0-s->bias_ratchet_strain/capacity,0.0);
    return p->bias_ratchet_rate*stress_level*bias_factor*saturation;
}

RIVA_HD static inline double riva_bias_ratchet_increment(
    const riva_parameters_t *p,const riva_state_t *old,const riva_state_t *trial,
    riva_tensor_t *direction)
{
    *direction=riva_zero();
    const double rate=riva_bias_ratchet_rate_factor(p,old);
    const double dl=riva_max(trial->lambda_total-old->lambda_total,0.0);
    const double norm=riva_norm(old->cyclic_direction);
    const double projection=riva_ddot(old->static_bias_tensor,old->cyclic_direction);
    if (rate<=0.0 || dl<=0.0 || norm<=1.0e-12) return 0.0;
    const double sign=projection>0.0?1.0:(projection<0.0?-1.0:0.0);
    *direction=riva_scale(old->cyclic_direction,sign/norm);
    const double capacity=riva_bias_ratchet_capacity(p,old);
    return riva_min(rate*dl,riva_max(capacity-old->bias_ratchet_strain,0.0));
}

RIVA_HD static inline riva_state_t riva_without_reversible_bias(
    const riva_parameters_t *p,const riva_material_parameters_t *material,
    const riva_state_t *state)
{
    if (fabs(state->bias_reversible_volume)<=1.0e-16) return *state;
    riva_state_t mechanical=*state;
    mechanical.eps_v_total=state->physical_eps_v_total;
    mechanical.eps_v_confining=mechanical.eps_v_total-
        mechanical.eps_v_irreversible-mechanical.eps_v_reversible;
    riva_tensor_t dev=riva_dev(state->stress);
    int at_floor=0;
    const double pressure=riva_pressure_from_confining(
        p,material,mechanical.eps_v_confining,
        mechanical.pressure_anchor,&at_floor); (void)at_floor;
    if (p->p_residual>0.0 || p->geostatic_admission_enabled) {
        dev=riva_regularize_deviator(
            p,material,pressure,mechanical.void_ratio,state->alpha,
            dev,&mechanical.alpha);
    } else mechanical.alpha=riva_scale(dev,1.0/pressure);
    mechanical.stress=riva_sub(
        dev,riva_iso(riva_physical_pressure(p,pressure)));
    return mechanical;
}

RIVA_HD static inline riva_state_t riva_forward_euler(
    const riva_parameters_t *p,const riva_material_parameters_t *material,
    const riva_state_t *old,riva_tensor_t deps,int force_reversal,int allow_legacy)
{
    const double physical_deps_v=riva_trace(deps);
    const riva_state_t mechanical=riva_without_reversible_bias(p,material,old);
    const riva_state_t trial=riva_backbone_forward_euler(
        p,material,&mechanical,deps,force_reversal,allow_legacy);
    riva_tensor_t ratchet_direction;
    const double ratchet=riva_bias_ratchet_increment(p,old,&trial,&ratchet_direction);
    riva_state_t provisional=trial;
    if (ratchet>0.0) {
        const riva_tensor_t effective_deps=riva_sub(deps,riva_scale(ratchet_direction,ratchet));
        provisional=riva_backbone_forward_euler(
            p,material,&mechanical,effective_deps,
            force_reversal,allow_legacy);
    }
    if (provisional.amplitude_reversals>old->amplitude_reversals) {
        const riva_tensor_t reversal_dev=riva_dev(old->stress);
        const riva_tensor_t excursion=riva_sub(reversal_dev,old->last_reversal_deviator);
        const double norm=riva_norm(excursion);
        if (norm>1.0e-12) provisional.cyclic_direction=riva_scale(excursion,1.0/norm);
    }
    const double target=riva_bias_reversible_volume_target(p,&provisional);
    riva_state_t state=provisional;
    state.bias_ratchet_strain=old->bias_ratchet_strain+ratchet;
    state.physical_eps_v_total=old->physical_eps_v_total+physical_deps_v;
    state.bias_reversible_volume=target;
    const double effective_total=state.physical_eps_v_total+target;
    const int changed=fabs(effective_total-state.eps_v_total)>1.0e-16;
    state.eps_v_total=effective_total;
    state.eps_v_confining=state.eps_v_total-state.eps_v_irreversible-
        state.eps_v_reversible;
    if (changed) {
        riva_tensor_t dev=riva_dev(state.stress); int at_floor=0;
        const double pressure=riva_pressure_from_confining(
            p,material,state.eps_v_confining,
            state.pressure_anchor,&at_floor);
        if (p->p_residual>0.0 || p->geostatic_admission_enabled) {
            dev=riva_regularize_deviator(
                p,material,pressure,state.void_ratio,old->alpha,
                dev,&state.alpha);
        } else state.alpha=riva_scale(dev,1.0/pressure);
        state.stress=riva_sub(
            dev,riva_iso(riva_physical_pressure(p,pressure)));
        state.pressure_floor_hits += at_floor?1:0;
    }
    state.void_ratio=old->void_ratio+(1.0+old->void_ratio)*physical_deps_v;
    state.D_re *= riva_bias_reversible_factor(p,&state);
    state.D=state.D_ir+state.D_re;
    return state;
}

RIVA_HD static inline riva_tensor_t riva_host_direction(const riva_parameters_t *p,
    riva_tensor_t deps,int *valid)
{
    const riva_tensor_t dev=riva_dev(deps); const double norm=riva_norm(dev);
    *valid=norm>p->reversal_strain_deadband;
    return *valid?riva_scale(dev,1.0/norm):riva_zero();
}

RIVA_HD static inline int riva_host_reversal(const riva_parameters_t *p,
    const riva_state_t *s,riva_tensor_t direction,int direction_valid)
{
    if (!direction_valid || !s->cyclic_phase_active) return 0;
    const double previous_norm=riva_norm(s->last_host_deviatoric_strain_direction);
    if (previous_norm<=1.0e-14) return 0;
    const double cosine=riva_ddot(direction,s->last_host_deviatoric_strain_direction)/previous_norm;
    if (cosine>p->reversal_direction_cosine) return 0;
    const double excursion=riva_norm(riva_sub(riva_dev(s->stress),s->last_reversal_deviator));
    return excursion>=p->reversal_stress_deadband_ratio*
        riva_max(s->pressure_anchor,riva_cone_pressure_floor(p));
}

/* reversal_override: -1 = detect from the committed state and this increment
 * (the original behavior, bit-identical); 0/1 = force the reversal decision.
 * The override exists so an implicit host can LATCH the decision made on the
 * first trial of an equilibrium increment and reuse it for later Newton
 * iterations of the same increment: re-detecting per trial makes the response
 * branch flip between iterations, and near a high static shear bias the two
 * branches differ enough that Newton limit-cycles instead of converging. */
RIVA_HD static inline int riva_update_material_ex(const riva_parameters_t *p,
    const riva_material_parameters_t *material,riva_tensor_t deps,
    int32_t fixed_substeps,riva_state_t *state,riva_tensor_t *stress_new,
    riva_update_info_t *info,int32_t reversal_override)
{
    if (!state || !state->initialized || fixed_substeps<1 ||
        !riva_material_parameters_valid(p,material) ||
        !riva_finite_tensor(deps)) return 0;
    /* A rate-independent zero increment must not advance compatibility,
     * reversal, or diagnostic history merely because the host took a step. */
    if (riva_ddot(deps,deps) == 0.0) {
        if (stress_new) *stress_new=state->stress;
        if (info) { info->accepted_substeps=0;
                    info->reversal_registered=0; }
        return 1;
    }
    int valid=0;
    const riva_tensor_t direction=p->objective_reversal_enabled?
        riva_host_direction(p,deps,&valid):riva_zero();
    const int reversal=p->objective_reversal_enabled?
        (reversal_override<0?riva_host_reversal(p,state,direction,valid):
         (reversal_override!=0)):0;
    riva_state_t current=*state;
    const riva_tensor_t sub=riva_scale(deps,1.0/(double)fixed_substeps);
    for (int32_t i=0;i<fixed_substeps;i++)
        current=riva_forward_euler(p,material,&current,sub,
            p->objective_reversal_enabled && reversal && i==0,
            !p->objective_reversal_enabled);
    if (p->objective_reversal_enabled && valid)
        current.last_host_deviatoric_strain_direction=direction;
    if (!riva_finite_tensor(current.stress) ||
        !isfinite(current.lambda_total) || !isfinite(current.void_ratio)) return 0;
    *state=current; if (stress_new) *stress_new=current.stress;
    if (info) { info->accepted_substeps=fixed_substeps;
                info->reversal_registered=reversal; }
    return 1;
}

RIVA_HD static inline int riva_update_material(const riva_parameters_t *p,
    const riva_material_parameters_t *material,riva_tensor_t deps,
    int32_t fixed_substeps,riva_state_t *state,riva_tensor_t *stress_new,
    riva_update_info_t *info)
{
    return riva_update_material_ex(p,material,deps,fixed_substeps,state,
                                   stress_new,info,-1);
}

RIVA_HD static inline int riva_update(const riva_parameters_t *p,riva_tensor_t deps,
    int32_t fixed_substeps,riva_state_t *state,riva_tensor_t *stress_new,
    riva_update_info_t *info)
{
    const riva_material_parameters_t material=
        riva_reference_material_parameters(p);
    return riva_update_material(p,&material,deps,fixed_substeps,state,
                               stress_new,info);
}

RIVA_HD static inline double riva_compatibility_residual(const riva_state_t *s)
{ return s->eps_v_total-s->eps_v_confining-s->eps_v_irreversible-s->eps_v_reversible; }

RIVA_HD static inline int riva_state_values(
    const riva_state_t *state,double values[RIVA_STATE_VALUE_COUNT])
{
    int i=0;
#define RIVA_STATE_TENSOR(t) do { \
    values[i++]=(t).xx; values[i++]=(t).yy; values[i++]=(t).zz; \
    values[i++]=(t).xy; values[i++]=(t).yz; values[i++]=(t).xz; \
} while (0)
    RIVA_STATE_TENSOR(state->stress); RIVA_STATE_TENSOR(state->alpha);
    RIVA_STATE_TENSOR(state->alpha0); RIVA_STATE_TENSOR(state->alpha01);
    RIVA_STATE_TENSOR(state->n); RIVA_STATE_TENSOR(state->fabric);
    values[i++]=state->D; values[i++]=state->beta;
    values[i++]=state->lambda_total; values[i++]=state->ep_eq_since_reversal;
    values[i++]=state->void_ratio; values[i++]=(double)state->reversals;
    values[i++]=(double)state->pressure_floor_hits;
    values[i++]=(double)state->denominator_floor_hits;
    values[i++]=(double)state->beta_fallbacks;
    values[i++]=state->eps_v_total; values[i++]=state->eps_v_confining;
    values[i++]=state->eps_v_irreversible; values[i++]=state->eps_v_reversible;
    values[i++]=state->pressure_anchor; values[i++]=state->D_ir;
    values[i++]=state->D_re; RIVA_STATE_TENSOR(state->last_reversal_deviator);
    values[i++]=state->cyclic_amplitude; values[i++]=state->amplitude_factor;
    values[i++]=(double)state->amplitude_reversals;
    values[i++]=state->initial_relative_state;
    values[i++]=state->state_contraction_factor;
    values[i++]=state->effective_knee_ratio;
    RIVA_STATE_TENSOR(state->geostatic_deviator);
    RIVA_STATE_TENSOR(state->static_bias_tensor);
    RIVA_STATE_TENSOR(state->cyclic_direction);
    values[i++]=state->static_bias_index;
    values[i++]=(double)state->cyclic_phase_active;
    values[i++]=state->bias_ratchet_strain;
    values[i++]=state->physical_eps_v_total;
    values[i++]=state->bias_reversible_volume;
    RIVA_STATE_TENSOR(state->last_host_deviatoric_strain_direction);
    values[i++]=state->ep_half_last;
    values[i++]=state->eta_m;
#undef RIVA_STATE_TENSOR
    return i;
}

#undef RIVA_HD
#endif /* OPENSEES_RIVA_SAND_KERNEL_H */
