# RIVA-Sand intermediate-bias research material

This directory exercises the separately named OpenSees research material
`RIVASandIntermediateBiasResearch`. It is the native implementation of the
constitutive successor developed on
`research/rivasand-intermediate-bias-plastic-gate`.
It does not replace or modify the frozen production `RIVASand` material.

This branch appends a completed-half-cycle plastic-activity memory to the
research state. It gates only the inherited reversible bias wave, and dynamic
activation does not preload a fictitious reversal count. This is kernel
restart revision 3 and is intentionally incompatible with revision-2 restart
files.

## Command

```tcl
nDMaterial RIVASandIntermediateBiasResearch tag Dr M kd h m zeta \
    eMax eMin Q R nG \
    <-rho value> <-nSub integer> <-stressScale value> \
    <-pMin value> <-tangentPMin value> <-pResidual value> \
    <-geostaticAdmission> <-reversalLatch> <-stage 0|1|2> \
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
```

The standalone replay in `tests/RIVASandIntermediateBiasResearchNativeReplay.cpp`
compares the allocation-free native kernel against the private six-history
handoff oracle without using Python at runtime.

The restart test enables `-reversalLatch` and verifies that both an element
material copy and a material reconstructed through `sendSelf`/`recvSelf`
retain the setting. The standalone state-contract test additionally checks
that automatic detection is unchanged and that a forced reversal is applied
exactly once across fixed constitutive substeps.

The standalone state-contract test can be built without OpenSees libraries:

```sh
c++ -std=c++17 -O2 \
  EXAMPLES/RIVASandIntermediateBiasResearch/tests/RIVASandIntermediateBiasResearchKernelStateTest.cpp \
  -o RIVASandIntermediateBiasResearchKernelStateTest
./RIVASandIntermediateBiasResearchKernelStateTest
```
