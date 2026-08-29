# RIVA-Sand intermediate-bias research material

This directory exercises the separately named OpenSees research material
`RIVASandIntermediateBiasResearch`. It is the native implementation of the
constitutive successor developed on `research/rivasand-intermediate-bias-flow`.
It does not replace or modify the frozen production `RIVASand` material.

## Command

```tcl
nDMaterial RIVASandIntermediateBiasResearch tag Dr M kd h m zeta \
    eMax eMin Q R nG \
    <-rho value> <-nSub integer> <-stressScale value> \
    <-pMin value> <-tangentPMin value> <-pResidual value> \
    <-geostaticAdmission> <-stage 0|1> \
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
substeps per host strain increment. As with production RIVA-Sand, use stage 0
for gravity and activate all element copies with:

```tcl
updateMaterialStage -material $matTag -stage 1
```

Add `-geostaticAdmission` only when the gravity workflow intentionally admits
a compressive geostatic state outside the cyclic bounding surface.

## Tests

After building OpenSees, run:

```sh
./build-riva-ib/OpenSees \
  EXAMPLES/RIVASandIntermediateBiasResearch/tests/RIVASandIntermediateBiasResearch_stage_activation.tcl

./build-riva-ib/OpenSees \
  EXAMPLES/RIVASandIntermediateBiasResearch/RIVASandIntermediateBiasResearch_material_point.tcl
```

The standalone replay in `tests/RIVASandIntermediateBiasResearchNativeReplay.cpp`
compares the allocation-free native kernel against the private six-history
handoff oracle without using Python at runtime.
