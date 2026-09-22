# RIVASAND02 OpenSees material

`RIVASAND02` is the OpenSees name of the intermediate-bias successor to
`RIVASand`, previously named `RIVASandIntermediateBiasResearch`. The `02`
identifies this successor; it does not refer to the historical PJ-Liq V2 model.
The original `RIVASand` material remains available separately.

The opt-in `-reversalType` selects the earlier, loading-branch-reference,
or host-increment UMAT-type reversal rule. Types 2 and 3 are research modes.
The separately named built-in `RIVASAND02BranchReversalResearch` remains
available and is equivalent to `RIVASAND02 -reversalType 2`.
See [column numerical corrections and validation](COLUMN_NUMERICS_VALIDATION.md)
for its command, restrictions, regression tests, retained SSPbrickUP failure
propagation, and the removal of the committed-stiffness damping change.

The old command RIVASandIntermediateBiasResearch is no longer registered.
The current command requires twelve material values: insert dimensionless
G0 immediately after Dr (reference: 483.48301127222084). Pre-G0 inputs and databases
are incompatible; rerun initialization. The 139-value kernel history remains
revision 4, but its adapter now serializes G0 plus an input revision marker
(183 values total). Existing stage behavior and BiasVolume options are unchanged.
Custom eMax/eMin also refresh the material-specific reference density,
matching the corrected Hercules implementation.
The read-only material response `referenceRelativeDensity` reports that
derived calibration anchor; it survives material copies and database restart.

See [G0 input migration and verification](G0_INPUT_VALIDATION.md) for the
completed checks and the regression command for both material generations.

This rename preserves the constitutive equations, calibration, and accepted
histories from `research/rivasand-field-bias-volume`; it does not expand their
validation scope. The reported limitations of the field correction (mode 2)
still apply.

This branch appends a completed-half-cycle plastic-activity memory to the
research state. It gates only the inherited reversible bias wave, and dynamic
activation does not preload a fictitious reversal count. The optional field
bias-volume correction described below makes this kernel restart revision 4;
it is intentionally incompatible with revision-3 restart files.

## Command

```tcl
nDMaterial RIVASAND02 tag Dr G0 M kd h m zeta \
    eMax eMin Q R nG \
    <-rho value> <-nSub integer> <-stressScale value> \
    <-pMin value> <-tangentPMin value> <-TanType 0|1> <-pResidual value> \
    <-geostaticAdmission> <-reversalType 1|2|3> <-reversalLatch> \
    <-BiasVolume 0|1|2> \
    <-stage 0|1|2> \
    <-initialStress sxx syy szz sxy syz sxz>
```

The public inputs and units follow [the RIVA-Sand user guide](../RIVASand/RIVASand_USER_GUIDE.md), including its G0 calibration and migration instructions. The additional calibrated
phase-transformation, loose-flow, mapping/backstress, and intermediate-bias
controls are frozen inside this research kernel; they are not independent
OpenSees inputs.

The Ottawa F65 reference values used by the verification cases are:

```tcl
set Dr 0.662
nDMaterial RIVASAND02 8001 \
    $Dr 483.48301127222084 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 \
    -nSub 1 -stressScale 1.0
```

Use `-stressScale 1.0` when the OpenSees stress unit is kPa and
`-stressScale 1000.0` when it is Pa. `-nSub` selects fixed constitutive
substeps per host strain increment.

`-reversalType` is an integer selected at material creation:

| Value | Reversal rule |
|---|---|
| `1` (default) | Earlier model: compare the trial strain-increment direction with the preceding committed increment direction |
| `2` | Reversal-reference research model: compare against accumulated loading-branch strain |
| `3` | UMAT-type research model: evaluate the stress-ratio predictor rule once for the full trial increment |

All three make one decision before constitutive substepping and apply any
reversal reset in the first substep. With the latch disabled, each Newton
trial recomputes that decision from the committed material state. This does
not freeze the decision across Newton iterations. Type 3 preserves the
tested host-increment variant's strain-based initialization until cyclic
activation, followed by its literal UMAT-type rule; it adds no switching
tolerance or convergence guard.

Existing `RIVASAND02` inputs retain type 1. The older
`RIVASAND02BranchReversalResearch` command retains type 2 and accepts only
`-reversalType 2` if the option is supplied explicitly. Identical repeated
selections are allowed; conflicting, missing, fractional, and out-of-range
values are rejected. The selection cannot be changed with `updateParameter`.

Material copies retain the selection. Query it using
`eleResponse $eleTag reversalType` for `SSPbrickUP`, or
`eleResponse $eleTag material 1 reversalType` for `bbarBrick`.
Type 1 retains its existing database format and latch support. Types 2 and 3
reject `-reversalLatch` and disable checkpoint/channel transfer; they must
be rerun from initialization. Their calibration, mesh/timestep objectivity,
and coupled-analysis convergence remain research limitations. Selecting
a reversal rule does not select or modify element damping.

`-TanType 0` (default) preserves the elastic tangent. `-TanType 1` selects
the safeguarded continuum elastoplastic backbone tangent, including the
mapping/backstress branch. It does not change stress integration or any
BiasVolume mode, and it is not a fully consistent algorithmic tangent of
the research overlays/substeps. Stage 0 and the initial tangent remain
elastic. Use a nonsymmetric-capable system solver. The selector and
committed tangent-activity flag survive copying/restart in adapter revision
2; the 139-value kernel state remains revision 4. The reader also accepts
the preceding G0-aware adapter-revision-1 checkpoints as elastic mode 0.
See [usage and safeguards](../RIVASand/RIVASand_USER_GUIDE.md#tangent-and-convergence-considerations)
and [validation](../RIVASand/CONTINUUM_TANGENT_VALIDATION.md).

`-BiasVolume` is a dimensionless integer selector for the three existing
biased-volume behaviors. Omitting it is exactly equivalent to `-BiasVolume 0`.

| Value | Behavior | Equivalent previous input |
|---|---|---|
| `0` (default) | Original research response: inherited reversible biased-volume wave and mean component enabled, field correction off | Neither old flag |
| `1` | Disable the inherited reversible biased-volume wave and mean component; keep the phase-transformation and irreversible channels | `-noBiasVolume` |
| `2` | Retain the inherited wave and enable the selective field correction to its mean component | `-fieldBiasVolume` |

The selector chooses a mode at material creation; it is not a scale factor
and cannot be changed with `updateParameter` during an analysis. The old
flags remain supported as aliases. Repeating the same selection is allowed,
including `-BiasVolume 1 -noBiasVolume`; requesting different modes is an
error regardless of argument order. Values outside 0--2, fractional values,
and missing values are rejected.

The selected mode survives material copying and database `save`/`restore`
through the existing revision-4 configuration flags. Query it for a
`bbarBrick` integration point with:

```tcl
eleResponse $eleTag material 1 BiasVolume
```

For `SSPbrickUP`, use `eleResponse $eleTag BiasVolume`. The existing
`noBiasVolume` and `fieldBiasVolume` boolean responses are also retained.

`-reversalLatch` is available only with reversal type 1. It is an opt-in OpenSees iteration stabilizer for dynamic
research analyses. It makes one host-level reversal decision on the first
accepted material evaluation of a load step and reuses that decision during
subsequent Newton trial evaluations. The transient decision is cleared on
commit, revert, and restart; the enabled setting is preserved by material
copying and `sendSelf`/`recvSelf`. The option is disabled by default, so the
accepted constitutive histories remain unchanged. Because the decision is
tied to the first accepted Newton trial, timestep and iteration-path
objectivity must still be demonstrated for a complete boundary-value problem
before treating it as a production default.

`-BiasVolume 2` (alias `-fieldBiasVolume`) is an opt-in constitutive research
correction for sloping-ground boundary-value analyses. It continuously removes only the contractive
part of the mean component of the inherited reversible biased-volume target
as the state enters the low-static-bias, high-confinement field window. A
positive, dilative mean shift is never reduced. The correction is bounded to
the inherited contractive mean magnitude, so it cannot reverse that
component. It leaves the oscillatory pressure wave, phase transformation,
shear ratchet, stress mapping, and plastic flow unchanged. The smooth windows
are internally fixed at static-
bias index 0.20--0.28 and pressure-anchor ratio 1.00--1.25 relative to the
40-kPa mean-transition pressure. A second smooth limiter fades the correction
out as the current effective-pressure ratio falls from 0.35 to 0.10. This
limiter was intended to reduce the correction near the cone apex, but it
does not establish stability of the coupled boundary-value analysis. Its
activity grows monotonically with plastic multiplier on a fixed 0.0001 scale. Activity
is committed only after an accepted host increment and is held fixed through
all constitutive substeps and repeated Newton trials from that committed
state. This is a research option, not part of the production `RIVASand`
calibration.

Mode 2 remains experimental: the UCD_6 comparison reported in
[PR #9](https://github.com/zwyll/OpenSees_RIVA/pull/9) found slower convergence
and worse pore-pressure agreement than mode 1. The unified selector preserves
those existing behaviors; it does not correct that reported limitation.

`-BiasVolume 1` (alias `-noBiasVolume`) is a separate opt-in diagnostic option
carried forward from the Phase-3 field study. It disables the entire inherited reversible
biased-volume target, including both its oscillatory pressure wave and mean
shift. It does not disable the phase-transformation or irreversible
dilatancy/contraction channels. The option remains off by default and is not
a recalibrated constitutive replacement for the biased-volume law. OpenSees
rejects a command that combines `-noBiasVolume` with `-fieldBiasVolume`, since
there is no mean biased-volume component left for the field correction to
modify.

For a controlled research comparison using the field correction, select:

```tcl
lappend matCmd -geostaticAdmission -reversalLatch -BiasVolume 2
```

For a diagnostic rerun that completely disables the inherited biased-volume
target, use instead:

```tcl
lappend matCmd -geostaticAdmission -reversalLatch -BiasVolume 1
```

The `-fieldBiasVolume` option retains the calibrated oscillatory pressure
wave; `-noBiasVolume` disables the whole reversible biased-volume target.
They cannot be combined. Select `-BiasVolume 0` (or omit the selector) for
the original research response. `-reversalLatch` remains an independent,
opt-in option.

Without `-geostaticAdmission`, the conventional two-stage sequence is
unchanged: use stage 0 for gravity and activate the nonlinear cyclic material
with:

```tcl
updateMaterialStage -material $matTag -stage 1
```

Add `-geostaticAdmission` only when the gravity workflow intentionally admits
a compressive geostatic state outside the cyclic bounding surface. In that
workflow, three stages keep geostatic equilibration separate from the
calibrated cyclic research mechanisms:

```tcl
# Stage 0: elastic gravity loading
updateMaterialStage -material $matTag -stage 1
# Stage 1: nonlinear geostatic admission and drained re-equilibration
# ...complete the gravity settle, permeability switch, and quiet hold...
updateMaterialStage -material $matTag -stage 2
# Stage 2: mapping/backstress and phase-transformation mechanisms active
# ...begin the dynamic analysis...
```

Stage 1 retains the converged skeleton stresses and uses the non-expansive
geostatic-admission rule until each over-bound stress point re-enters the
ordinary cone. Stage 2 is a stress-preserving state transition; it must be
issued before dynamic loading when `-geostaticAdmission` is enabled. This
staging change introduces no new calibrated material parameter.

## Tests

After building OpenSees, run:

```sh
./build-riva-ib/OpenSees \
  EXAMPLES/RIVASAND02/tests/RIVASAND02_stage_activation.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASAND02/RIVASAND02_material_point.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASAND02/tests/RIVASAND02_restart.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASAND02/tests/RIVASAND02_no_bias_volume.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASAND02/tests/RIVASAND02_bias_volume_modes.tcl
```

The standalone replay in `tests/RIVASAND02NativeReplay.cpp`
compares the allocation-free native kernel against the private six-history
handoff oracle without using Python at runtime.

The restart test enables `-reversalLatch` and `-fieldBiasVolume` and verifies
that both an element material copy and a material reconstructed through
`sendSelf`/`recvSelf` retain both settings. The standalone state-contract test
additionally checks that automatic detection is unchanged, that a forced
reversal is applied exactly once across fixed constitutive substeps, and that
the field correction removes only the bounded mean component.

The no-bias-volume test verifies that the switch is off by default and that
the enabled setting survives both the element material copy and
`sendSelf`/`recvSelf`.

The mode-selector test compares cyclic stress and complete state histories
for modes 0, 1, and 2 against their previous inputs, checks restart continuation
for every mode, and rejects malformed or conflicting mode selections.

The standalone state-contract test can be built without OpenSees libraries:

```sh
c++ -std=c++17 -O2 \
  EXAMPLES/RIVASAND02/tests/RIVASAND02KernelStateTest.cpp \
  -o RIVASAND02KernelStateTest
./RIVASAND02KernelStateTest
```

The state-contract test also checks once-per-increment reversal scheduling
at nSub = 1, 4, 16, and 40 for all three reversal types. The column regression
runner includes `reversal_types.tcl`, which checks selector parsing, default
and alias equivalence, material copies, full-state rollback, replacement of
speculative trials, tangent-choice independence of prescribed stress/state
histories, and the research-mode restrictions. See
[the column validation instructions](COLUMN_NUMERICS_VALIDATION.md) to run
it with a built OpenSees executable and the diagnostic probe. These checks
verify the implementation contracts; they do not establish convergence of
the sloping-ground analysis.

The selector integration was also compared with the saved earlier and
branch-reference implementations at commit `e97c12fbc`, and with the saved
literal host-increment UMAT research build. Each rule covered 24 prescribed
rotating/changing-pressure cases: Pa/kPa, direct cyclic activation or
geostatic-to-dynamic activation, all three BiasVolume modes, and both tangent
types. Across 6,240 accepted states per rule, types 1 and 2 matched exactly.
Type 3 matched all stresses and the first 138 state entries exactly; its
remaining activity entry differed by at most `1.11e-16`. No full-duration
slope convergence or new parameter calibration is claimed by this comparison.
