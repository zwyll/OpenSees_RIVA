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
1.11e-16. No new full-duration free-field or CDSS fit is claimed by this port.

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
