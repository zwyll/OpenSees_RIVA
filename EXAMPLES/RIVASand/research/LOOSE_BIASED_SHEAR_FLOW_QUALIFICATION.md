# Loose-biased shear-flow qualification

Date: 2026-08-27

Research branch: `research/rivasand-loose-biased-shear-flow`

Constitutive checkpoint before qualification: `44b4f5430`

## Decision

The loose-biased shear-flow successor is **not qualified to replace production
RIVA-Sand**.  It remains a useful research checkpoint: it corrects the
`4666_loose_b015` calibration history and improves several broad shape and
pore-pressure measures relative to frozen production.  However, it fails all
18 untouched hard gates, terminates before four prescribed comparison cycles,
and extrapolates poorly from static-bias ratio 0.15 to 0.25--0.375.

The production `research/rivasand-geostatic` branch remains unchanged.  A
native successor port was deliberately not created after this failed
constitutive gate.  The native production kernel was compiled, replayed, and
benchmarked to establish the baseline that a future accepted successor must
beat or closely approach.

## Qualification scope

The qualification contains:

- the original eight-case PRJ-3484/PRJ-4666 audit at 32 points with four local
  substeps and 64 points with two local substeps;
- a frozen, untouched 18-case holdout matrix covering both projects, four
  density bands, 40 and 100 kPa confinement, static-bias ratios from zero to
  0.375, and CSR from 0.15 to 0.50;
- direct comparison against frozen production and the signed-PT parent;
- 32-by-4 versus 64-by-2 host/local resolution checks;
- disabled and out-of-window equivalence, restart, finiteness, compatibility,
  state-copy, and local-substep tests; and
- compiled native production golden replay and constitutive runtime tests.

The holdout manifest validates every source file, scalar test condition, and
experimental 7.5%-strain crossing before running a prediction.  No holdout
parameter was recalibrated.

## Original eight-case audit

Seven of eight cycles-to-triggering predictions meet the factor-1.35 gate at
32-by-4; all eight meet it at 64-by-2.  None passes the complete hard gate,
because loop-shape or pore-pressure-history requirements remain outside their
limits.

| Case | Experiment N | 32-by-4 N | 64-by-2 N |
|---|---:|---:|---:|
| 3484_u015 | 32.164 | runout | 26.156 |
| 3484_u021 | 4.180 | 4.531 | 3.531 |
| 3484_b025 | 66.227 | 65.250 | 67.203 |
| 3484_b030 | 3.219 | 3.188 | 3.188 |
| 4666_loose_b015 | 7.094 | 7.188 | 7.188 |
| 4666_dense_u038 | 11.258 | 11.531 | 11.656 |
| 4666_dense_b025 | 36.234 | 33.250 | 34.250 |
| 4666_dense_b0375 | 29.188 | 29.219 | 30.219 |

The target `4666_loose_b015` result is stable from 32 through 128 host points:
`N = 7.188`, `7.188`, and `7.180` at 32-by-4, 64-by-2, and 128-by-2.  The
known 16-by-8 coarse limit gives `N = 6.250` and remains outside the practical
objectivity range.

## Untouched holdout audit

| Case | Experiment N | Production | Signed-PT parent | Successor 32-by-4 | Successor 64-by-2 |
|---|---:|---:|---:|---:|---:|
| 3484 Dr55, unbiased, CSR 0.15 | 10.719 | 11.031 | 11.031 | 11.031 | 9.031 |
| 3484 Dr71, unbiased, CSR 0.17 | 17.227 | runout | runout | runout | 20.578 |
| 3484 Dr76, unbiased, CSR 0.21 | 10.156 | runout | runout | runout | 15.547 |
| 3484 100 kPa, unbiased, CSR 0.17 | 7.156 | 8.531 | 8.531 | 8.531 | 7.203 |
| 3484 bias 0.25, CSR 0.25 | 14.180 | 5.094 | 11.219 | 11.219 | 11.234 |
| 3484 bias 0.375, CSR 0.35 | 14.203 | 6.188 | runout | runout | runout |
| 4666 Dr50, unbiased, CSR 0.15 | 7.711 | 6.031 | 6.031 | 6.031 | 5.031 |
| 4666 Dr60, unbiased, CSR 0.15 | 12.695 | 10.531 | 10.531 | 10.531 | 9.031 |
| 4666 Dr90, unbiased, CSR 0.35 | 18.219 | 14.531 | 10.531 | 10.531 | 10.547 |
| 4666 Dr50, bias 0.15, CSR 0.18 | 20.109 | 3.094 | 25.250 | 28.219 | 28.219 |
| 4666 Dr50, bias 0.15, CSR 0.21 | 10.117 | 0.688 | 7.063 | 7.188 | 7.203 |
| 4666 Dr60, bias 0.15, CSR 0.22 | 12.133 | 2.156 | 5.094 | 5.031 | 5.109 |
| 4666 Dr50, bias 0.075, CSR 0.17 | 10.125 | 6.125 | 15.031 | runout | runout |
| 4666 Dr50, bias 0.25, CSR 0.30 | 9.055 | 1.156 | 5.250 | 2.219 | 2.188 |
| 4666 Dr60, bias 0.25, CSR 0.33 | 7.141 | failed | 2.188 | 1.188 | 2.203 |
| 4666 Dr90, bias 0.25, CSR 0.45 | 23.242 | 5.188 | 4.250 | 4.250 | 4.250 |
| 4666 Dr50, bias 0.375, CSR 0.35 | 13.188 | 5.219 | 19.219 | 1.219 | 1.172 |
| 4666 Dr90, bias 0.375, CSR 0.50 | 29.188 | 17.219 | 21.219 | 21.219 | 21.234 |

Summary at 32-by-4:

| Model | Joint score | Comparison cycles available | N within factor 1.35 | Complete hard gates |
|---|---:|---:|---:|---:|
| Frozen production | 6.654 | 10/18 | 5/18 | 0/18 |
| Signed-PT parent | 3.711 | 16/18 | 6/18 | 0/18 |
| Shear-flow successor | 4.481 | 14/18 | 5/18 | 0/18 |

The successor improves normalized full-history `r_u` error from 1.917 for
production and 1.560 for the signed-PT parent to 1.294.  It also improves raw
loop-phase error to 0.571 from 0.787 and 1.057, respectively.  Those gains do
not compensate for the four unavailable loops or the high-bias CSR--N
failures.

## Resolution and state tests

All five inherited and successor state/regression programs pass.  The direct
successor strain-history check reports a 6.351% maximum stress difference
between two and four local substeps.  Restart is exact, compatibility remains
closed, all states remain finite, and disabling the mechanism or leaving its
declared density/bias window recovers the parent exactly.

The broader holdout is not timestep-objective enough for acceptance.  Seven
cases change by more than 10% or switch between a finite trigger and runout
between 32-by-4 and 64-by-2.  The largest finite discrepancy occurs for the
Dr60, bias-0.25, CSR-0.33 case (`N = 1.188` versus `2.203`).  Several of these
differences are inherited outside the new correction window, but a production
successor must still resolve them as a complete model.

## Native production baseline

The native C++ production kernel passes the five frozen golden histories plus
reference/custom equivalence, state-copy isolation, restart, zero-increment,
and finite-state tests.  The benchmark used GCC 15.1.0, `-O3 -DNDEBUG`, arm64,
seven measured samples, and 192,000 host updates per sample.

| Fixed substeps | Median time (s) | ns/host update | ns/local substep |
|---:|---:|---:|---:|
| 1 | 0.203 | 1,056 | 1,056 |
| 2 | 0.386 | 2,010 | 1,005 |
| 4 | 0.750 | 3,908 | 977 |
| 8 | 1.482 | 7,719 | 965 |

The approximately constant cost per local substep confirms linear scaling of
the frozen native kernel.  The Python research oracle previously measured the
successor at 1.52 times its signed-PT parent, but that number is not a native
runtime ratio.  Reporting it as one would be misleading.  A true native
successor benchmark is deferred until the constitutive holdout gate passes;
otherwise the port would create a second implementation of a rejected model.

## Required next development

The next mechanism should address two coupled failures before any native port:

1. replace the loose correction's monotonic extrapolation above bias 0.15,
   which causes premature triggering at bias 0.25--0.375, with a physically
   smooth bias-dependent hardening/ratchet balance; and
2. remove the inherited host-resolution sensitivity seen in the unbiased and
   intermediate-density holdouts without changing the accurate early-cycle
   target loop.

The untouched matrix must remain frozen while that work proceeds.  Passing
only `4666_loose_b015` is no longer an acceptable promotion criterion.

## Reproduction

```text
python RIVASandLooseBiasedShearFlowHoldoutAudit.py \
  --output results/qualification/holdout_successor_32x4 \
  --model successor --points 32 --substeps 4 --full

python RIVASandLooseBiasedShearFlowHoldoutAudit.py \
  --output results/qualification/holdout_successor_64x2 \
  --model successor --points 64 --substeps 2 --full

g++ -std=c++17 -O3 -DNDEBUG RIVASandNativeKernelBenchmark.cpp \
  -o RIVASandNativeKernelBenchmark
./RIVASandNativeKernelBenchmark
```
