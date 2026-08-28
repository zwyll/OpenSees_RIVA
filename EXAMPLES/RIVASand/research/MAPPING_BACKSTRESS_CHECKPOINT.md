# Directional mapping-backstress checkpoint

Date: 2026-08-28

Branch: `research/rivasand-mapping-backstress`

Starting checkpoint: `ef8d809d8` on
`research/rivasand-loose-biased-shear-flow`

## Decision

Retain this implementation as a **numerical research checkpoint**, but do not
promote or port it to production OpenSees or Hercules yet.  It demonstrates
that a directional backstress inside the bounding-surface mapping plus a
fixed stress corrector can remove the severe local-substep bifurcation of the
rejected scalar overlay.  It also corrects cycles to the strain criterion for
the three low/intermediate-density high-bias development cases.

It does not yet pass the complete constitutive gate: those three numerical
loops remain too narrow and too dissipative, their pressure-wave histories
remain inaccurate, and one case retains material host-timestep sensitivity.
Frozen production RIVA-Sand and the qualified predecessor remain unchanged.

## Isolated activation

The mechanism is active only for low/intermediate-density sand above the
preserved alpha=0.15 static-bias state.  Disabled mode, alpha=0.15, the Ottawa
3484 b025 case, and the dense cases reproduce every predecessor state field
exactly.  The original eight-case 32-by-4 results are therefore unchanged:

| Case | Experiment N | Predecessor and checkpoint N |
|---|---:|---:|
| 3484 u015 | 32.164 | runout |
| 3484 u021 | 4.180 | 4.531 |
| 3484 b025 | 66.227 | 65.250 |
| 3484 b030 | 3.219 | 3.188 |
| 4666 loose b015 | 7.094 | 7.188 |
| 4666 dense u038 | 11.258 | 11.531 |
| 4666 dense b025 | 36.234 | 33.250 |
| 4666 dense b0375 | 29.188 | 29.219 |

The full 18-case holdout audit likewise reproduces the predecessor outside
the three declared active cases.

## Active high-bias results

At the reference 32 host points per cycle and four fixed local substeps:

| Case | Experiment N | Predecessor N | Mapping checkpoint N |
|---|---:|---:|---:|
| 4666 Dr50, alpha 0.25, CSR 0.30 | 9.055 | 2.219 | 9.156 |
| 4666 Dr60, alpha 0.25, CSR 0.33 | 7.141 | 1.188 | 7.125 |
| 4666 Dr50, alpha 0.375, CSR 0.35 | 13.188 | 1.219 | 12.156 |

The corresponding prescribed-cycle loop metrics expose the unresolved shape
error:

| Case | Range exp/model (%) | Area exp/model | Damping exp/model | Ru RMSE |
|---|---:|---:|---:|---:|
| Dr50, alpha 0.25 | 1.405 / 0.792 | 17.207 / 14.317 | 0.328 / 0.481 | 0.143 |
| Dr60, alpha 0.25 | 2.575 / 0.841 | 27.740 / 15.897 | 0.262 / 0.458 | 0.248 |
| Dr50, alpha 0.375 | 1.252 / 0.639 | 16.191 / 14.767 | 0.297 / 0.527 | 0.094 |

Thus the triggering-cycle fit cannot be treated as a successful joint
calibration.  The stress--strain loops still fail the frozen hard gates.

## Numerical gates

The direct state test passes:

- exact disabled and preserved-window equivalence;
- finite response and exact state-copy/restart replay;
- reversible--irreversible volume compatibility;
- outer-cone stress correction;
- a monotone local shear-stress map; and
- 2/4 and 4/8 fixed-substep stress differences of 2.768% and 1.176% in the
  non-liquefied local integration probe.

Stress-controlled resolution results are:

| Case | 32x4 N | 64x2 N | 128x2 N | 128x4 N |
|---|---:|---:|---:|---:|
| Dr50, alpha 0.25 | 9.156 | 9.141 | 9.148 | 9.156 |
| Dr60, alpha 0.25 | 7.125 | 6.188 | 6.180 | 6.188 |
| Dr50, alpha 0.375 | 12.156 | 12.156 | 11.094 | 11.086 |

The near identity of 128x2 and 128x4 confirms that the fixed local corrector
removed local-substep dependence.  The remaining 32-to-128 difference is a
host-step dependence in the inherited phase-volume update, not a local
stress-root bifurcation.  It is still too large for production acceptance in
the Dr60 case.

## Cost

A five-sample Python material-point timing with 1,000 host updates and four
fixed substeps gave median times of 3.228 s for the predecessor and 6.201 s
for this checkpoint, a ratio of 1.92.  This is an interpreter-level research
ratio, not a native C++ or GPU prediction.  No native implementation should be
created until the constitutive gate passes.

## Next constitutive requirement

The stress corrector and translated tensor mapping should be retained.  The
next change must address loop width and energy together through a recoverable
directional memory surface (or an equivalent branch-elastic domain) inside
the mapping.  It must reduce plastic area while increasing recoverable strain
range; changing contraction, scalar modulus, or ratchet capacity again would
only tune N without fixing damping.  The inherited host-level phase-volume
integration must also be made semigroup-consistent before a production port.

## Reproduction

```text
python RIVASandMappingBackstressTest.py

python RIVASandMappingBackstressAudit.py \
  --output results/mapping_backstress/original8_32x4 \
  --points 32 --substeps 4 --full

python RIVASandMappingBackstressHoldoutAudit.py \
  --output results/mapping_backstress/holdout18_32x4 \
  --points 32 --substeps 4 --full
```
