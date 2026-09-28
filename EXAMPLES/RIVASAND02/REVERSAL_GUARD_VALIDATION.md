# Opt-in type-3 reversal eligibility guard

Branch `research/rivasand02-reversal-guard` starts from
`d88a1123d90c065bedc9f7f30f0ac652281bee19` on
`research/rivasand02-column-numerics`. It exposes the previously isolated
small-loading guard as a per-material option. The default remains off.
This is an unaccepted RIVASAND02 research successor event rule.

## Usage

Append these options to an existing complete material definition:

```tcl
-reversalType 3 -reversalGuard 0.0001
```

`-reversalGuard 0` or omission preserves the literal type-3 criterion.
Positive values require type 3; zero is allowed with every type. The
threshold is dimensionless and uses the material's admitted pressure scale,
so it does not need rescaling between Pa and kPa models. It controls the
minimum deviatoric-stress excursion accumulated since the last reversal.
The inherited tiny-strain cutoff is enabled with the positive guard.

The guard makes one eligibility decision before local material substeps.
Every Newton trial recomputes from committed history. It introduces no
Newton latch, new history slots, damping change or calibration change.
Type 3 continues to reject `-reversalLatch` and checkpoint/channel transfer.
Type 1 retains the existing restart format. Material copies and resets
preserve the guard setting; different materials may use different values
within the same domain. Query it with `eleResponse $eleTag reversalGuard`
for SSPbrickUP, or `eleResponse $eleTag material 1 reversalGuard` for bbarBrick.

## Verification

Run the synthetic adapter suite with a built executable and the diagnostic
probe described by `tests/column_numerics/run.py`:

```sh
python EXAMPLES/RIVASAND02/tests/column_numerics/run.py \
  --opensees /path/to/OpenSees --probe /path/to/columnprobe.so \
  --output /path/to/guard-checks
c++ -std=c++17 -O2 EXAMPLES/RIVASAND02/tests/RIVASAND02ReversalGuardTest.cpp \
  -o /tmp/rivasand02-guard-test
/tmp/rivasand02-guard-test
```

The guard checks cover invalid inputs, material copying, reset, coexistence
of guarded and unguarded materials, suppression of tiny-cycle events,
retention of finite-cycle reversals, rollback and replacement of speculative
trials, and at most one event per host increment at nSub 1/4/16/40. The
existing tests continue to check both tangent options, the other reversal
rules, stage activation, bias-volume options, G0 input, type-1 restart,
SSPbrickUP material-failure propagation and trial-stiffness damping.

Native build comparisons use prescribed paths across Pa/kPa, geostatic
activation, BiasVolume 0/1/2 and TanType 0/1. The guard-off paths are compared
with the parent executable; guard-on paths are compared with the archived
prototype. These implementation checks do not establish coupled convergence.

The native branch build passed both standalone kernel tests, all five
column scripts and all five existing adapter scripts. Guard-off histories
were byte-identical to the parent for each reversal type: 24 cases and
6,240 accepted states per type. At guard 0.0001, all stresses, the first
138 state entries and tangent-status values matched the archived prototype
exactly in another 24 cases; the last activity entry differed by at most
1.11e-16. These port checks did not establish a new calibration or CDSS fit.

## Native free-field initialization check (2026-09-28)

A subsequent fixed-parameter three-layer free-field comparison uses the
current type-1 calibrated coefficients, `-nSub 40`, `-TanType 0`, no latch,
and NewmarkExplicit 0.5 / Linear. Gravity and the final ten settling steps
remain elastic. Nonlinear activation is followed by five seconds with no
earthquake excitation, using the dynamic permeability and damping settings.
The earthquake is admitted only after the equilibrium check passes, with
the committed material state and pressure integral preserved.

At both 0.00125 and 0.000625 s, types 1 and 2 and guarded type 3 have zero
new reversals during the hold. All 139 material-state entries are identical
across those passing variants at each timestep, both before and after the
hold. Maximum relative shear-stiffness drift is below 4e-10. Unguarded type 3
instead accumulates 2,916 / 4,539 reversals and retains only 15% of its initial
shear stiffness. The guard-off executable reproduces the parent hold
byte-for-byte. Elastic settling alone therefore does not remove the
unguarded type-3 startup defect in this case.

The initialization screen checks all nodes and elements at every step:
displacement 1e-7 m, velocity 1e-7 m/s, acceleration 1e-5 m/s2, pressure and
stress change 1e-4 kPa, strain change 1e-8, zero new reversals, and relative
shear-stiffness change 1e-4. These limits screen stationary initialization;
they are not a coupled-response convergence criterion.

Type 2 and guarded type 3 each complete the subsequent 70-second measured
earthquake at both timesteps, with no failed steps or material rejections.
Spectra use the native acceleration samples. With the finer run as reference:

| Rule | Largest sensor PGA difference | Largest sensor Sa log-RMS difference | Largest pressure-peak difference | Largest normalized pressure-history RMS difference |
|---|---:|---:|---:|---:|
| Type 2 | 11.452% | 0.108386 | 0.277% | 0.432% |
| Type 3, guard 0.0001 | 7.564% | 0.065729 | 0.720% | 0.290% |

Both fail the all-sensor screen (5% for PGA and pressure, ln(1.05) for
spectral log RMS), with the largest acceleration discrepancy at 8.7 m depth.
Experimental fit remains similar; no parameters were refitted. Linear
explicit stepping does not provide a nonlinear residual-convergence test.
The guard is useful for this initialization procedure, but this evidence
does not establish earthquake timestep independence or a new accepted
calibration. Raw histories and the detailed research report remain local.

## Limits

The prior diagnostic study found fewer near-zero reversal events and
improved stationary initialization for selected steps. Full-duration
implicit convergence remained unresolved; the finest explicit timestep
comparison narrowly failed its sensor-PGA screen. Some calibrated biased
CDSS paths also changed before their stress-control failure, and measured
free-field fit did not improve. Consequently `0.0001` is a tested research
setting, not a recommended default or a newly accepted calibration.

Failure propagation and the restored trial-stiffness damping are inherited
unchanged from the parent branch. The earlier committed-stiffness damping
patch is not reintroduced. Hercules and original RIVASand are unchanged.
