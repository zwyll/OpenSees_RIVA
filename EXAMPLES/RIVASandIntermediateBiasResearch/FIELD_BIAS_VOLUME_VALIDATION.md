# Field-bias mean correction audit

## Status

`-fieldBiasVolume` is an opt-in research mechanism. It is not part of
production RIVA-Sand and is not qualified for general field analysis. The
implemented form is deliberately narrow: it can remove only a contractive
mean shift at low static bias and at an admitted mean pressure above the
40-kPa sign-transition pressure. It never suppresses a positive dilative mean
shift.

Its activity is a bounded committed state variable. Plastic activity advances
it once after a successful host increment; the value is fixed through all
constitutive substeps and repeated Newton trials from that committed state.
This avoids an iteration-dependent material switch.

## Preserved material-point histories

With the option disabled, all legacy columns in the six accepted native
histories are byte-exact relative to the parent
`research/rivasand-intermediate-bias-reversal-latch` checkpoint:

- PRJ-3484 zero bias
- PRJ-3484 intermediate bias 0.25
- PRJ-3484 intermediate bias 0.30
- PRJ-4666 loose bias 0.15
- PRJ-4666 dense zero bias
- PRJ-4666 loose high bias

These cases do not calibrate the new high-confinement window. They demonstrate
only that the opt-in mechanism does not change the accepted default response.

## RPI-A boundary-value audit

The RPI-A target file gives P1 `r_u,max = 0.88169` and end-of-motion
`r_u = 0.69368`. The unmodified 25 by 10 parent model already gives 0.89433
and 0.69166, respectively. Its maximum admitted element mean pressure is
28.24 kPa, below the 40-kPa mean-shift transition; therefore the final
high-confinement correction is exactly inactive for this mesh.

Two broader trial forms were rejected:

| trial | mesh | P1 `r_u,max` | P1 end `r_u` | recovered cuts | peak surface displacement | wall time |
|---|---:|---:|---:|---:|---:|---:|
| parent checkpoint | 25 by 10 | 0.894 | 0.692 | 32 | 0.177 m | 776 s |
| cancel the mean term below 40 kPa | 25 by 10 | 0.820 | 0.467 | 44 | 0.196 m | 771 s |
| parent checkpoint | 10 by 4 | 0.753 | 0.209 | 0 | 0.050 m | 58 s |
| cancel the mean term below 40 kPa | 10 by 4 | 0.901 | 0.168 | 1 | 0.044 m | 57 s |
| add one mean-term magnitude below 40 kPa | 10 by 4 | 0.839 | 0.112 | 0 | 0.044 m | 57 s |

Both low-pressure corrections reduce the retained P1 pore pressure. The
reason is constitutive rather than numerical: below 40 kPa the inherited mean
shift is positive and dilative. Removing it, or amplifying it without jointly
changing the pressure-wave phase, does not improve the coupled field response.

## Remaining qualification

No local RPI-B input motion is present in this workspace, so an RPI-B run was
not performed. A future active calibration of the final high-confinement
window requires an untouched low-bias case whose admitted mean pressure
actually exceeds 40 kPa. Until such data are available, keep
`-fieldBiasVolume` off in production analyses.
