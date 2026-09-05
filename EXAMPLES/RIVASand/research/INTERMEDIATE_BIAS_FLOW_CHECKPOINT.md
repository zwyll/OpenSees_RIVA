# Intermediate-density biased-flow and compatible-volume research checkpoint

## Scope and ancestry

This private research successor starts from mapping/backstress commit
`411afa994`. It does not change production RIVA-Sand. The low-density
mapping/backstress update, restrained phase-transformation volume law,
accumulation control, loose-biased correction, and all production source under
`SRC/` are retained.

The first target was the branch-to-branch strain imbalance in the two PRJ-3484
intermediate-density biased tests. A substep study established that the open
loops were neither local nonconvergence nor false reversal detection. The
second target was the inherited one-sided pressure wave in `3484_b030` and the
false 7.5% DA trigger in the 100-kPa, CSR=0.20 runout.

## Added mechanism

The successor adds a reversal-anchored directional shear-compliance term to
the modulus used inside the existing local stress update. It is:

- centered on the `Dr=66.3%` density band;
- zero at zero bias and outside compact low- and high-bias windows;
- continuous through its density, bias, amplitude, branch-progress, and
  phase-zone coordinates;
- directionally balanced through the fixed admitted static-bias direction;
- inactive in the qualified low-density mapping/backstress window and in the
  dense Ottawa histories; and
- state-free apart from two immutable density--bias gate values cached at
  cyclic activation.

The high-bias branch now also:

- evaluates its amplitude gate on a consistent full-cycle scale at the first
  reversal;
- activates smoothly at the lower stabilized amplitude of the CSR=0.20
  runout instead of applying a one-step first-reversal impulse;
- shifts the reversible PT reference from the one-sided admitted static
  potential toward the center of the cyclic potential range;
- admits the compatible PT volume continuously before the first formal
  reversal from the stress excursion away from the admitted static state;
- ramps that PT activity over completed half-cycles; and
- slows only the high-bias reversible-volume relaxation.

No `r_u` multiplier or direct pressure correction is used. The reversible and
irreversible volumes are committed once per host increment and effective
pressure is rebuilt from compatible volumetric strain. The bounding surface,
reversal detector, and low-bias shear calibration are unchanged.

## Calibrated controls

| Control | Value |
|---|---:|
| intermediate density onset / peak / cutoff | 0.60 / 0.663 / 0.76 |
| low-bias onset / full / cutoff onset / cutoff | 0.30 / 0.50 / 0.60 / 0.66 |
| low-bias compliance peak / bell gain | 1.00 / 0.50 |
| low-bias directional balance | 0.00 |
| high-bias onset / full / cutoff onset / cutoff | 0.62 / 0.67 / 0.72 / 0.78 |
| high-bias amplitude onset / full | 0.40 / 0.45 |
| high-bias compliance peak / bell gain | 3.00 / 1.50 |
| high-bias directional balance | 0.12 |
| minimum intermediate branch multiplier | 0.20 |
| first-reversal amplitude scale | 0.40 |
| high-bias PT anchor fraction | 0.90 |
| high-bias PT relaxation multiplier | 0.25 |
| high-bias PT activation length | 6 reversals |
| pre-reversal PT activity scale | 1.50 |

## Targeted 32-point / 4-substep results

| Case | Quantity | Previous checkpoint | Successor | Experiment |
|---|---|---:|---:|---:|
| `3484_b025` | cycles to criterion | 65.250 | 66.219 | 66.227 |
|  | cycle-11 center (%) | 3.389 | 3.266 | 3.337 |
|  | range ratio | 0.654 | 0.981 | 1.000 |
|  | area ratio | 0.916 | 1.133 | 1.000 |
|  | damping ratio | 0.297 | 0.245 | 0.212 |
|  | phase-normalized loop RMSE | 0.164 | 0.111 | 0.000 |
| `3484_b030` | cycles to criterion | 3.188 | 3.219 | 3.219 |
|  | cycle-2 center (%) | 6.018 | 5.443 | 5.361 |
|  | range ratio | 0.566 | 1.015 | 1.000 |
|  | area ratio | 0.669 | 0.865 | 1.000 |
|  | damping ratio | 0.374 | 0.270 | 0.316 |
|  | phase-normalized loop RMSE | 0.363 | 0.045 | 0.000 |

The compatible-volume correction is exactly inactive in `3484_b025`; all 83
pre-existing history fields remain bit-exact. In `3484_b030`, it preserves
the calibrated shear loop and changes the pressure metrics as follows:

| `3484_b030` quantity | Shear-flow checkpoint | Coupled-volume successor | Experiment |
|---|---:|---:|---:|
| cycles to criterion | 3.219 | 3.219 | 3.219 |
| cycle-2 center (%) | 5.443 | 5.438 | 5.361 |
| strain-range ratio | 1.015 | 1.015 | 1.000 |
| loop-area ratio | 0.865 | 0.864 | 1.000 |
| phase-normalized loop RMSE | 0.0449 | 0.0448 | 0.0000 |
| full-history `r_u` RMSE | 0.1764 | 0.1182 | 0.0000 |
| cycle-2 `r_u` mean error | 0.0409 | 0.0018 | 0.0000 |
| cycle-2 `r_u` amplitude ratio | 0.919 | 0.970 | 1.000 |
| cycle-2 `r_u` phase error | 21.0 deg | 10.9 deg | 0.0 deg |

`3484_b030` now passes the strict combined shear, pressure-history, and cycle
gate. Its early first-cycle pressure transient remains the largest residual.

## Broader checks

- The six original cases outside the intermediate gate retain their prior
  32x4 cycle predictions exactly.
- The three low-density, high-bias PRJ-4666 holdouts retain their prior cycle
  predictions exactly: 9.156, 7.125, and 12.156 cycles.
- The affected 40-kPa, alpha=0.25 CSR--N log RMSE changes only from 0.288 to
  0.284. The existing overly steep local segment near CSR 0.27 remains.
- At 100 kPa and alpha=0.30, CSR=0.15, 0.20, and 0.25 all remain runouts for
  their complete experimental durations; the prior false CSR=0.20 trigger at
  13.25 cycles is eliminated. CSR=0.30 remains at 3.219 cycles.
- The seven original histories outside the new high-bias volumetric gate are
  bit-exact over all 83 common recorded fields. The three selected low-density
  PRJ-4666 event counts remain 9.156, 7.125, and 12.156 cycles.
- With the new mechanism disabled, every inherited mapping/backstress state
  field is bit-exact after a multibranch cyclic strain path.
- Focused tests verify non-overlapping bias windows and exact inactivity for
  zero-bias and dense states.
- `3484_b030` reaches 3.219 cycles with two, four, and eight fixed constitutive
  substeps. Its `r_u` RMSE is 0.117, 0.118, and 0.119, respectively.

## Computational cost and decision

Before port preparation, a direct Python strain-driven benchmark with four
fixed substeps showed measurable wrapper and active-mechanism overhead. A
response-preserving optimization pass then removed the following redundant
work:

- scalar smoothstep coordinates no longer dispatch through NumPy array
  clipping;
- initial and reference relative density are evaluated once instead of on
  every constitutive query;
- zero intermediate and loose-flow gates return before evaluating branch
  progress; and
- the inherited base, phase-transformation, and loose cyclic-flow factors
  share one pressure activity rather than evaluating it three times.

For ten strain-driven cycles at 32 host points and four fixed substeps, the
median of five Python timings changed from 3888 to 2184 microseconds per host
update in the `3484_b025` window, from 5711 to 3336 in the `3484_b030` window,
and from 2215 to 1207 on a zero-bias inactive path. These are reductions of
44%, 42%, and 46%, respectively. They are Python-oracle measurements, not
native OpenSees or Hercules predictions; scalar NumPy overhead will not exist
in the native port.

The optimization is response-exact in the available regression record. All
eight primary histories are bit-for-bit identical over 87 recorded fields,
the three high-bias holdout histories are bit-for-bit identical, and all 11
affected CSR jobs retain identical non-runtime records.

### Native-port evaluation cache

The flattened Hercules kernel should preserve these reductions and go one
step further by using an update-local scalar cache. Evaluate the following
once at initialization or cyclic activation: initial relative density,
reference relative density, density-only gates, admitted static-bias gates,
and the fixed static-bias direction. Evaluate the following once per trial
state/substep and share them among modulus, hardening, flow, and compatible
volume laws: stress invariants, pressure ratio, bounding/dilatancy surfaces,
projected bias, branch progress, transformation zone, and cyclic-flow
activity. Recorder-only diagnostics must remain outside the constitutive
inner loop.

The ratchet-corrected second backbone evaluation, fixed constitutive
substeps, host-level phase-volume update, and final compatible-pressure
rebuild are not redundant. They change the accepted state and must remain in
the native implementation.

The correction resolves the two stated defects in the research driver, but it
is not yet a production replacement. The inherited 40-kPa alpha=0.25 CSR--N
local-slope error, the `3484_b025` shear hard-gate failure, the residual first-
cycle `3484_b030` pressure transient, and native OpenSees/Hercules runtime and
restart qualification remain before any port or production rename.
