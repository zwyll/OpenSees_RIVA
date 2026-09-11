# Optional continuum backbone tangent — 2026-09-11

Both `RIVASand` and `RIVASAND02` accept `-TanType 0|1`. Omitted/0 preserves
the elastic default; 1 requests the safeguarded continuum elastoplastic
backbone tangent. Required inputs, stress integration, calibration, kernel
histories, and the 02 BiasVolume/reversalLatch settings are unchanged.
The implementation is on `research/rivasand-continuum-tangent`, based on
`feature/rivasand` at `2b8d374ba`.

This is an approximate continuum operator with frozen auxiliary overlays,
not the consistent algorithmic derivative of the full update. No finite-
difference stress replays are performed during analysis. See the
[user guide](RIVASand_USER_GUIDE.md#tangent-and-convergence-considerations)
for nonsymmetric solvers, fallbacks, U–P ownership, damping, and limitations.

## Completed checks

- Rebuilt the complete local OpenSees executable with GNU C++ 15 on macOS.
- Six-column smooth-backbone finite-difference rate check: contraction,
  dilation, and zero dilatancy in Pa/kPa. Maximum bulk-normalized error at
  the smallest perturbation was below `1.7e-7`. This test intentionally
  does **not** claim agreement with the complete finite-substep Jacobian.
- All five frozen original native histories and all six successor native
  histories pass; the largest successor normalized discrepancy against the
  existing oracle was approximately `3.6e-11`. No golden files or kernel
  equations were changed.
- Original and successor native state/admission tests; existing Tcl stage,
  restart, G0, no-bias-volume and three-mode BiasVolume tests pass.
- New prescribed-strain tests exercise both materials, 1/4/16 substeps,
  zero/intermediate/high bias and the 02 mapping branch, and all BiasVolume
  modes. Stresses and complete kernel state histories are **exactly equal**
  for omitted, 0, and 1 tangent selections. Read-only tangent queries do not
  advance history.
- Default adapter histories also match the separately retained pre-change
  executable byte-for-byte, not just the new executable's explicit mode 0.
- Tests cover invalid/missing/fractional/conflicting options, per-element
  copies, Pa/kPa tangents, pressure-floor fallback, nonlinear restart, and
  gravity/admission/stage transitions with continuum selected.
- New G0-aware adapter revision 2 preserves vector sizes (135 / 183) and
  kernel layouts. It stores the selector and committed plastic-activity
  flag in adapter configuration bits. The new reader successfully resumes
  actual databases written by the preceding revision-1 executable as mode 0,
  with identical continuation. Old executables reject revision-2 databases;
  pre-G0 databases remain incompatible.

## Small boundary-value probes

The paired tests use identical loads, fixed host steps, four constitutive
substeps, Newton, BandGeneral, and a displacement-increment tolerance of
`1e-10`. There is no Rayleigh damping. U–P runs use `brickUP`, finite water
compressibility, and three seconds of loading at `dt=0.01 s`; they are not
LEAP mesh or full-motion tests. The biased U–P case has initial shear of
10 kPa at vertical effective stress 40 kPa.

| Material | Probe | Completed, both modes | Newton iterations 0 → 1 | Median wall seconds 0 → 1 |
|---|---|---|---|---|
| RIVASand | Monotonic static ramp | 128/128 | 535 → 371 | 0.171 → 0.124 |
| RIVASAND02 | Monotonic static ramp | 128/128 | 535 → 371 | 0.197 → 0.140 |
| RIVASand | Unbiased U–P dynamics | 300/300 | 1482 → 906 | 0.570 → 0.363 |
| RIVASAND02 | Unbiased U–P dynamics | 300/300 | 1482 → 906 | 0.689 → 0.440 |
| RIVASand | Biased U–P dynamics | 300/300 | 1484 → 1278 | 0.623 → 0.568 |
| RIVASAND02 | Biased U–P dynamics | 300/300 | 1207 → 975 | 0.742 → 0.673 |

Times are medians of five local repetitions with the existing CMake build
(no Release build type), not optimized material-kernel throughput or a
prediction of full-mesh speedup. These short timings are indicative only.

**Retained failure:** cyclic force-controlled static probes at amplitudes
1 and 4 kPa stop on the first unloading step, after 8/128 accepted steps,
for both materials and both tangent modes. Mode 1 does not cure this existing
reversal/global-solve limitation. The benchmark prints these failed cases
deliberately; successful process exit does not mean every probe converged.

Mode 1 is therefore an opt-in solver aid, not a new production-qualification
claim. Full slice/LEAP meshes, stronger motion, timestep sensitivity, and
stiffness-proportional damping still require problem-specific validation.
Changing Newton iterations can also change a first-trial reversal latch
decision; prescribed-strain parity alone does not guarantee identical
finite-tolerance boundary-value trajectories.

## Reproduce

From the repository, compile/run the standalone rate test:

```sh
c++ -std=c++17 -O2 EXAMPLES/RIVASand/tests/RIVASandContinuumRateTest.cpp -o /tmp/RIVASandContinuumRateTest
/tmp/RIVASandContinuumRateTest
```

Run these with the rebuilt OpenSees executable from a scratch working
directory, using absolute script paths:

```text
EXAMPLES/RIVASand/tests/RIVASand_continuum.tcl
EXAMPLES/RIVASand/tests/RIVASand_continuum_stages.tcl
EXAMPLES/RIVASand/tests/RIVASand_tangent_benchmark.tcl
```

For independent default/restart comparisons, the continuum script accepts
an output filename and uses only old-compatible commands in that mode.
Compare output produced by the old and new executables. Run
`RIVASand_tangent_checkpoint.tcl make` with the old executable, then
`RIVASand_tangent_checkpoint.tcl check` with the new one in the same scratch
directory. These checks leave their comparison databases in that directory.
