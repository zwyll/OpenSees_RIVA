# RIVA-Sand intermediate-bias research material

This directory exercises the separately named OpenSees research material
`RIVASandIntermediateBiasResearch`. It is the native implementation of the
constitutive successor developed on
`research/rivasand-intermediate-bias-plastic-gate`.
It does not replace or modify the frozen production `RIVASand` material.

This branch appends a completed-half-cycle plastic-activity memory to the
research state. It gates only the inherited reversible bias wave, and dynamic
activation does not preload a fictitious reversal count. The optional field
bias-volume correction described below makes this kernel restart revision 4;
it is intentionally incompatible with revision-3 restart files.

## Command

```tcl
nDMaterial RIVASandIntermediateBiasResearch tag Dr M kd h m zeta \
    eMax eMin Q R nG \
    <-rho value> <-nSub integer> <-stressScale value> \
    <-pMin value> <-tangentPMin value> <-pResidual value> \
    <-geostaticAdmission> <-reversalLatch> <-fieldBiasVolume> <-noBiasVolume> \
    <-stage 0|1|2> \
    <-initialStress sxx syy szz sxy syz sxz>
```

The public inputs and units follow `RIVASand`. The additional calibrated
phase-transformation, loose-flow, mapping/backstress, and intermediate-bias
controls are frozen inside this research kernel; they are not independent
OpenSees inputs.

The Ottawa F65 reference values used by the verification cases are:

```tcl
set Dr 0.662
nDMaterial RIVASandIntermediateBiasResearch 8001 \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 \
    -nSub 1 -stressScale 1.0
```

Use `-stressScale 1.0` when the OpenSees stress unit is kPa and
`-stressScale 1000.0` when it is Pa. `-nSub` selects fixed constitutive
substeps per host strain increment.

`-reversalLatch` is an opt-in OpenSees iteration stabilizer for dynamic
research analyses. It makes one host-level reversal decision on the first
accepted material evaluation of a load step and reuses that decision during
subsequent Newton trial evaluations. The transient decision is cleared on
commit, revert, and restart; the enabled setting is preserved by material
copying and `sendSelf`/`recvSelf`. The option is disabled by default, so the
accepted constitutive histories remain unchanged. Because the decision is
tied to the first accepted Newton trial, timestep and iteration-path
objectivity must still be demonstrated for a complete boundary-value problem
before treating it as a production default.

`-fieldBiasVolume` is an opt-in constitutive research correction for sloping
ground boundary-value analyses. It continuously removes only the contractive
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
prevents the correction from driving an already-liquefied point past the cone
apex while allowing it to oppose subsequent pressure rebuilding. Its activity
grows monotonically with plastic multiplier on a fixed 0.0001 scale. Activity
is committed only after an accepted host increment and is held fixed through
all constitutive substeps and repeated Newton trials from that committed
state. This is a research option, not part of the production `RIVASand`
calibration.

`-noBiasVolume` is a separate opt-in diagnostic option carried forward from
the Phase-3 field study. It disables the entire inherited reversible
biased-volume target, including both its oscillatory pressure wave and mean
shift. It does not disable the phase-transformation or irreversible
dilatancy/contraction channels. The option remains off by default and is not
a recalibrated constitutive replacement for the biased-volume law. OpenSees
rejects a command that combines `-noBiasVolume` with `-fieldBiasVolume`, since
there is no mean biased-volume component left for the field correction to
modify.

For an explicit research rerun of the three-stage sloping-ground workflow, use:

```tcl
lappend matCmd -geostaticAdmission -reversalLatch -fieldBiasVolume
```

For a diagnostic rerun that completely disables the inherited biased-volume
target, use instead:

```tcl
lappend matCmd -geostaticAdmission -reversalLatch -noBiasVolume
```

The `-fieldBiasVolume` option retains the calibrated oscillatory pressure
wave; `-noBiasVolume` disables the whole reversible biased-volume target.
They cannot be combined.

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
  EXAMPLES/RIVASandIntermediateBiasResearch/tests/RIVASandIntermediateBiasResearch_stage_activation.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASandIntermediateBiasResearch/RIVASandIntermediateBiasResearch_material_point.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASandIntermediateBiasResearch/tests/RIVASandIntermediateBiasResearch_restart.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASandIntermediateBiasResearch/tests/RIVASandIntermediateBiasResearch_no_bias_volume.tcl
```

The standalone replay in `tests/RIVASandIntermediateBiasResearchNativeReplay.cpp`
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

The standalone state-contract test can be built without OpenSees libraries:

```sh
c++ -std=c++17 -O2 \
  EXAMPLES/RIVASandIntermediateBiasResearch/tests/RIVASandIntermediateBiasResearchKernelStateTest.cpp \
  -o RIVASandIntermediateBiasResearchKernelStateTest
./RIVASandIntermediateBiasResearchKernelStateTest
```
