# Restrained accumulation-control research checkpoint

## Status and purpose

This is a private research successor to the restrained phase-transformation
checkpoint. It is **not** production RIVA-Sand and does not alter the frozen
production kernel. Its one purpose is to delay the excessive late-cycle shear
ratchet in the intermediate-density biased PRJ-3484 tests while preserving the
already calibrated early-cycle loop and effective-stress path.

## Added state and update

At cyclic activation, the model stores the current accumulated plastic
multiplier as `phase_accumulation_lambda_anchor`. The cyclic plastic-work
memory is

`W = max(lambda_total - phase_accumulation_lambda_anchor, 0)`.

The activity is a smooth product of:

- a compact density window, active from `Dr=0.53` to `0.82` and equal to one
  at `Dr=0.663`;
- a smooth static-bias gate that vanishes at zero bias;
- a smooth plastic-work transition above an amplitude-dependent onset.

The onset decreases with established cyclic amplitude `A`:

`W0(A) = 0.030 (0.43 / max(A, 0.43))`.

The transition width is `0.006`. The hardening prefactor is multiplied by

`1 + 3 gA(A) activity`,

where `gA` decreases smoothly from `1` at `A=0.43` to `2/3` at `A=0.65`.
This supplies a gain of three for 3484_b025 and approximately two for
3484_b030.

The amplitude-weighted activity is committed exactly once at the end of each
host increment and is held constant during every internal constitutive
substep. This prevents the memory evolution count from changing when only the
number of fixed constitutive substeps changes.

## Calibrated outcomes

The 32-point-per-cycle, four-substep audit gives:

| Case | Experiment N | Prior PT N | Restrained N |
|---|---:|---:|---:|
| 3484_b025 | 66.2266 | 26.21875 | 65.25 |
| 3484_b030 | 3.219 | 1.0625 | 3.1875 |

For 3484_b025, shear strain, shear stress, vertical effective stress, and
`ru_vertical` are bit-for-bit identical to the prior PT checkpoint through
cycle 11. The first stress and effective-stress differences occur at cycle
11.21875, after the calibrated early loop. The cycle-11 loop center is
`3.3893%` versus `3.3369%` experimentally, and its area ratio remains `0.916`.

The six non-target histories—3484_u015, 3484_u021, 4666_loose_b015,
4666_dense_u038, 4666_dense_b025, and 4666_dense_b0375—are exactly unchanged
in every common output field. The full eight-case score improves from about
`2.239` for the prior PT checkpoint to `1.697` here. The strict legacy gates
still fail; those gates include known shape and pore-pressure deficiencies and
must not be represented as production qualification.

## Integration audit

Changing from 32 points/cycle with four substeps to 64 points/cycle with two
substeps changes cycles to criterion by:

| Case | 32x4 N | 64x2 N | Difference |
|---|---:|---:|---:|
| 3484_b025 | 65.25 | 67.203125 | +2.99% |
| 3484_b030 | 3.1875 | 3.1875 | 0.00% |
| 4666_dense_b025 | 33.25 | 34.25 | +3.01% |
| 4666_dense_b0375 | 29.21875 | 30.21875 | +3.42% |

A severe strain-driven four-versus-eight-substep test has a 13.85% maximum
stress difference. Part of this is inherited from the underlying restrained
PT checkpoint, so the present model remains a research checkpoint rather than
a production replacement.

For the same fixed 512-increment strain history with four constitutive
substeps, a warm benchmark with three samples per model gave median update times of `1.156 s` for
the prior PT checkpoint and `1.183 s` for this model, or approximately `2.35%`
local-update overhead. A full-to-liquefaction b025 analysis naturally takes
longer because the corrected model now simulates about 65 rather than 26
cycles; that is additional physical history, not constitutive cost per step.

## Verification commands

From the OpenSees research worktree:

```text
python3 EXAMPLES/RIVASand/research/RIVASandAccumulationControlTest.py

MPLCONFIGDIR=/private/tmp/rivasand-mpl /Users/wzhang/anaconda3/bin/python \
  EXAMPLES/RIVASand/research/RIVASandAccumulationControlAudit.py \
  --output EXAMPLES/RIVASand/research/results/accumulation_control_host_committed_final_full \
  --points 32 --substeps 4 --full
```

## Decision boundary

This mechanism establishes that late-cycle restraint can fix the two
intermediate-density biased cycle counts without changing the early loop or
the other six audited histories. It should remain private research code until
its remaining loop-shape, pore-pressure-history, broader validation, and
integration-sensitivity limitations are resolved. It must not yet be ported
to production OpenSees or Hercules.
