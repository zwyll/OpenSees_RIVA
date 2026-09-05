# Loose-biased and dense-unbiased correction checkpoint

This private research successor starts from
`RIVASandAccumulationControlModel`. It is not the frozen production RIVA-Sand
model and must not be ported to OpenSees or Hercules without a separate
acceptance decision.

## Added mechanisms

1. **Loose biased stabilization.** For low-density sand under a committed
   geostatic static bias, the parent plastic modulus is multiplied by 1.40 and
   the irreversible contraction factor is multiplied by 0.50. Both changes are
   weighted by the same compact density/bias gate and vanish above
   `Dr = 0.56` or at zero static bias.
2. **Loose signed phase volume.** The branch-relative legacy pressure wave is
   replaced inside the same loose-biased gate. A continuous potential measured
   along the committed static-bias direction drives reversible volume with
   generic, wave, and mean activity scales of 0.20, 0.30, and 1.10. This removes
   the folded hourglass path without changing the calibrated triggering cycle.
3. **Dense zero-bias phase volume.** A fixed cyclic direction is captured after
   the first physical reversal. A smooth stress-ratio-magnitude potential then
   drives symmetric reversible phase volume on both half-cycles. The activity
   scale is 0.70 and the potential center is 0.80. The gate vanishes below
   `Dr = 0.82` and for committed static-bias index above 0.08.

Both gates use the committed geostatic `static_bias_index`. They do not use the
projected bias, which is undefined before the first cyclic direction is known.

## Eight-test result

The primary audit uses 32 target points per cycle and four fixed constitutive
substeps.

| Case | Experimental N | Numerical N | Result |
|---|---:|---:|---|
| 3484_u015 | 32.1641 | runout | inherited limitation |
| 3484_u021 | 4.17969 | 4.53125 | inherited |
| 3484_b025 | 66.2266 | 65.25 | inherited |
| 3484_b030 | 3.219 | 3.1875 | inherited |
| 4666_loose_b015 | 7.09375 | 7.125 | corrected |
| 4666_dense_u038 | 11.2578 | 11.53125 | corrected |
| 4666_dense_b025 | 36.2344 | 33.25 | inherited |
| 4666_dense_b0375 | 29.1875 | 29.21875 | inherited |

At 64 target points per cycle and two constitutive substeps, the two corrected
targets give `N = 7.109375` and `N = 11.65625`, respectively. The six out-of-window
32-by-4 histories reproduce the parent checkpoint to zero or roundoff-level
difference (`<= 1.6e-11` in the stored common fields).

## Remaining limitations

- The loose biased effective-stress path now has the experimental one-sided,
  smoothly curved topology instead of the prior folded hourglass. Its
  stress-strain loop remains too narrow and its center is too low.
- The dense unbiased effective-stress path now has the experimental broad
  two-lobed topology instead of the prior narrow twisted path, but the
  stress-strain loop remains wider than the experiment.
- The inherited 3484_u015 runout and 3484_u021 timestep sensitivity remain.
- The strict all-metric acceptance gate is not passed. This checkpoint is a
  targeted research improvement, not a production replacement.

## Reproduction

```text
python RIVASandLooseUnbiasedCorrectionTest.py
python RIVASandLooseUnbiasedCorrectionAudit.py \
  --output results/loose_unbiased_correction_final_full \
  --points 32 --substeps 4 --full
python plot_latest_all_datasets_comparison.py
```
