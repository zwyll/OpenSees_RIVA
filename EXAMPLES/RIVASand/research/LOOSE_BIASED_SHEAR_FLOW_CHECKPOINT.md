# Loose-biased shear-flow checkpoint

> **Qualification status (2026-08-27):** retained as a research checkpoint,
> but rejected as a production replacement after the untouched 18-case gate.
> See `LOOSE_BIASED_SHEAR_FLOW_QUALIFICATION.md`.

This private research successor starts from the signed phase-transformation
checkpoint in `RIVASandLooseUnbiasedCorrectionModel`. It is not the frozen
production RIVA-Sand model and must not be ported to OpenSees or Hercules
without a separate acceptance decision.

## Added mechanisms

The correction is multiplied by the existing compact low-density/static-bias
gate. It therefore activates for the PRJ-4666 loose biased test and vanishes
exactly for the other seven DSS histories.

1. **Late-branch shear compliance.** A smooth branch-progress coordinate
   reduces shear stiffness only near the end of each half-cycle. This enlarges
   the strain extrema while preserving the narrow mid-branch opening required
   for an S-shaped loop.
2. **Bounded directional ratchet.** Plastic work produces a finite strain shift
   along the committed geostatic-bias direction. The ratchet saturates at a
   tensorial strain of 0.0172 and does not reverse with the cyclic loading
   direction.
3. **Delayed plastic-work hardening.** The early-cycle hardening multiplier is
   1.40. A host-step-committed memory begins at accumulated plastic multiplier
   0.035 and approaches a late multiplier of 3.50 over a width of 0.010. This
   preserves the calibrated cycle-3 loop while preventing premature strain
   triggering.
4. **Signed PT retuning.** The already accepted continuous loose-sand pressure
   wave and mean activities are 0.45 and 1.40. The legacy branch-relative wave
   remains fully replaced, so the hourglass path does not return.

Additional calibrated shear controls are a cyclic modulus reduction of 0.85,
a loose shear-modulus scale of 0.75, a branch-compliance gain of 1.0 over branch
progress 0.60--0.95, and a directional-ratchet rate of 4.0.
The density/static-bias gate is evaluated once at cyclic activation and stored
in the material state; this does not change the response and avoids repeated
gate evaluation inside the local integration loop.

## Target result

The primary audit uses 32 target points per cycle and four fixed constitutive
substeps. Metrics are evaluated on cycle 3.

| Metric | Experiment | Signed-PT parent | Shear-flow successor |
|---|---:|---:|---:|
| Loop center (%) | 2.642 | 0.837 | 2.661 |
| Strain range (%) | 4.568 | 2.903 | 4.515 |
| Loop area (kPa-%) | 23.541 | 12.866 | 21.759 |
| Damping ratio | 0.180 | 0.153 | 0.167 |
| Mid-loop opening (%) | 0.309 | 0.399 | 0.405 |
| Mean `r_u` | 0.644 | 0.622 | 0.637 |
| `r_u` amplitude | 0.303 | 0.272 | 0.310 |
| `r_u` phase (degrees) | 160.7 | 163.3 | 161.6 |
| Triggering cycle | 7.094 | 7.125 | 7.188 |

The target-only joint audit score decreases from 0.511 to 0.411. Across all
eight cases, the joint score decreases from 1.208 to 1.192.

## Isolation and objectivity

- All seven non-target 32-by-4 and 64-by-2 histories are bit-exact to the
  signed-PT parent in every stored common field.
- The target triggering cycles are 7.1875 at 32-by-4, 7.1875 at 64-by-2,
  7.171875 at 64-by-4, and 7.1796875 at 128-by-2.
- At 64-by-2 the loop center, range, area, and damping ratio are 2.613%,
  4.499%, 21.385 kPa-%, and 0.166.
- The strain-driven regression gives a 6.351% maximum stress difference
  between two and four constitutive substeps and passes restart, compatibility,
  finiteness, disabled-equivalence, and out-of-window-equivalence checks.
- Three serial Python-oracle runs give median target-case runtimes of 3.329 s
  for the signed-PT parent and 5.062 s for this successor, a 1.52x ratio. A
  compiled native implementation should be benchmarked separately.

## Remaining limitations

- The numerical cycle-3 loop is still more angular than the experimental loop,
  and its mid-loop opening remains about 0.10 percentage point too large.
- Although the cycle-3 `r_u` mean, amplitude, and phase improve, the complete
  `r_u` history score is worse because the first numerical pressure peak is too
  large and later troughs remain too deep.
- A very coarse 16-point host cycle gives `N = 6.25`; the calibrated result is
  objective over the practical 32--128 point range, not at this coarse limit.
- The strict all-metric acceptance gate remains failed. This checkpoint is a
  targeted research improvement, not a production replacement.

## Reproduction

```text
python RIVASandLooseBiasedShearFlowTest.py
python RIVASandLooseBiasedShearFlowAudit.py \
  --output results/loose_biased_shear_flow_final_full \
  --points 32 --substeps 4 --full
python plot_latest_all_datasets_comparison.py
python plot_loose_biased_shear_flow_comparison.py
```
