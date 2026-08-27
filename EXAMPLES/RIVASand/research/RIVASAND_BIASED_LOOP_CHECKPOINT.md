# RIVA-Sand biased-loop research checkpoint

## Scope and ancestry

This private checkpoint starts exactly from
`origin/research/rivasand-geostatic` at commit `51d93d884`.  The OpenSees
`RIVASand` source under `SRC/material/nD/RIVASand/` is unchanged.  The new
Python files are an isolated research successor and must not be described or
ported as production RIVA-Sand without additional validation.

The rejected DBF/DMP/Iwan research worktree was not used as constitutive
source.  Its strict eight-case case matrix, stress servo, metrics, and gates
were reused only as a validation harness.

## Retained mechanism

The bounded calibration retained one compact change:

- a reversal-anchored Masing-type tangent-compliance multiplier;
- smoothly gated by dense-state weight and projected static bias;
- equal to one at a reversal and progressively softer across a half-cycle;
- a small branch balance limits excessive biased loop-center drift; and
- the existing cyclic shear-modulus reduction changes from 0.85 to 0.87.

It does **not** change the frozen hardening, dilatancy, compatibility,
pore-pressure, geostatic-admission, or bias-ratchet equations.  A proposed
signed stress-side gate was tested and removed because it improved one local
shape metric but produced unacceptable center drift.

Calibrated research controls are:

| Control | Value |
|---|---:|
| `cyclic_shear_modulus_reduction` | 0.87 |
| `branch_compliance_gain` | 4.5 |
| `branch_compliance_exponent` | 1.0 |
| `branch_compliance_bias_reference` | 0.5384061785684372 |
| `branch_compliance_bias_exponent` | 1.0 |
| `branch_compliance_minimum` | 0.10 |
| `branch_directional_balance` | 0.03 |
| `branch_balance_bias_exponent` | 0.20 |

## Fixed-cycle biased-loop results

The table compares production RIVA-Sand with the research checkpoint at each
experiment's declared comparison cycle.  Range and area are numerical divided
by experimental; a value of one is ideal.

| Case | Quantity | Production | Research | Experiment/target |
|---|---|---:|---:|---:|
| PRJ-3484, alpha=0.25, CSR=0.20 | loop center (%) | 3.375 | 3.269 | 3.337 |
|  | range ratio | 0.657 | 0.733 | 1.000 |
|  | area ratio | 0.911 | 0.890 | 1.000 |
|  | damping ratio | 0.294 | 0.257 | 0.212 |
| PRJ-4666 dense, alpha=0.25, CSR=0.38 | loop center (%) | 3.574 | 2.279 | 2.087 |
|  | range ratio | 0.210 | 0.775 | 1.000 |
|  | area ratio | 0.131 | 0.804 | 1.000 |
|  | damping ratio | 0.091 | 0.151 | 0.139 |
| PRJ-4666 dense, alpha=0.375, CSR=0.50 | loop center (%) | 4.809 | 2.648 | 2.568 |
|  | range ratio | 0.197 | 0.994 | 1.000 |
|  | area ratio | 0.126 | 1.209 | 1.000 |
|  | damping ratio | 0.081 | 0.154 | 0.127 |

The dense biased raw phase-normalized strain RMSE changes from 0.391 to 0.187
for alpha=0.25 and from 0.434 to 0.217 for alpha=0.375.  The corresponding
loop-center trajectory errors change from 0.363 to 0.161 and from 0.566 to
0.352.

## Full-cycle behavior

| Case | Experimental N | Production N | Research N |
|---|---:|---:|---:|
| PRJ-3484 alpha=0.25, CSR=0.20 | 66.23 | 27.16 | 27.22 |
| PRJ-4666 dense alpha=0.25, CSR=0.38 | 36.23 | 28.25 | 43.22 |
| PRJ-4666 dense alpha=0.375, CSR=0.50 | 29.19 | 24.22 | runout at 50.99 |

The alpha=0.25 dense cycle prediction improves.  The alpha=0.375 cycle result
becomes too resistant, so this checkpoint does not pass the CSR--N gate.
The pore-pressure histories are intentionally almost unchanged and retain the
production underprediction in both dense biased cases.

## Regression and numerical checks

- All five frozen production golden histories pass because production source
  is untouched.
- Kernel state/custom equivalence, trial isolation, restart, zero increment,
  residual-pressure translation, and stress-preserving geostatic admission
  pass.
- With the research mechanism disabled and the production value 0.85 restored,
  the research class is state-for-state identical to the standalone oracle.
- Research restart is exact in the Python state-copy test.
- A strain-driven 2-versus-4 substep check has 6.221% maximum stress difference.
- The 32-point/4-substep and 64-point/2-substep stress-controlled dense-loop
  ranges differ by less than 0.25%; loop areas differ by about 4%--6%.

## Decision

This is a useful stress--strain research checkpoint, not a production
replacement.  It fixes the most conspicuous dense biased narrow-ribbon and
center-drift errors with one inexpensive state-free modifier, but it does not
fully reproduce the measured S curvature and does not solve dense biased
pore-pressure or high-bias CSR--N behavior.  Do not port it to OpenSees or
Hercules production yet.
