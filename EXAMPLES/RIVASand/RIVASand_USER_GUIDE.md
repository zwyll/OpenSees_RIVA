# RIVA-Sand OpenSees user guide

## Purpose and scope

RIVA-Sand is a three-dimensional effective-stress `NDMaterial` for nonlinear
cyclic response and liquefaction analysis of saturated sand. The OpenSees
identifier is case-sensitive and is written `RIVASand`.

This document is intentionally limited to model use, input parameters,
initialization, output, and the calibrated Ottawa F65 parameter row. It does
not present the constitutive equations or other unpublished formulation
details.

## Material command

```tcl
nDMaterial RIVASand tag Dr M kd h m zeta eMax eMin Q R nG \
    <-rho value> <-nSub integer> <-stressScale value> \
    <-pMin value> <-tangentPMin value> <-stage 0|1> \
    <-initialStress sxx syy szz sxy syz sxz>
```

All eleven values from `Dr` through `nG` are required. OpenSees provides one
configurable `RIVASand` command; there is no separate `RIVASandCustom`
material.

## Ottawa F65 calibrated parameters

The following row is the frozen reference calibration used for the supplied
OpenSees implementation tests. All quantities in this table are
dimensionless.

| Input | Ottawa F65 value | Meaning | Practical effect |
|---|---:|---|---|
| `Dr` | 0.662962962962963 | Initial relative density, entered as a fraction | Defines the initial density state; this value corresponds to an initial void ratio of 0.601 |
| `M` | 1.25 | Base bounding stress ratio | Primarily controls the available shear-strength/stress-ratio range |
| `kd` | 1.125 | Base dilatancy stress ratio | Controls the stress-ratio level associated with the transition between contractive and dilative tendencies |
| `h` | 122.44207260468994 | Plastic hardening scale | Controls the magnitude of plastic hardening and therefore nonlinear stiffness |
| `m` | 0.945 | Bounding-surface mapping exponent | Controls how hardening changes as the stress state approaches its image point |
| `zeta` | 0.025 | Volumetric-coupling scale | Controls the magnitude of contractive and dilative volumetric response and the associated effective-stress evolution |
| `eMax` | 0.78 | Maximum/reference loose-state void ratio | Upper bound used to convert relative density to the initial void ratio |
| `eMin` | 0.51 | Minimum/reference dense-state void ratio | Lower bound used to convert relative density to the initial void ratio |
| `Q` | 10.0 | Pressure-dependent state-mapping coefficient | Works with `R` to define the reference density state as confinement changes |
| `R` | 1.5 | Pressure-dependent state-mapping coefficient | Works with `Q`; the two coefficients should normally be calibrated together |
| `nG` | 0.65 | Pressure exponent for elastic stiffness | Controls how small-strain stiffness changes with effective confinement |

`Dr` must be entered as a fraction between 0 and 1, not as a percentage. The
material converts it to its initial void ratio using `Dr`, `eMax`, and `eMin`.
For example, the reference inputs `Dr=0.662962962962963`, `eMax=0.78`, and
`eMin=0.51` give `e0=0.601`.

In OpenSees, the other ten required parameters do not change automatically
when `Dr` changes. For another relative density of Ottawa F65, the user may
retain the Ottawa values for `M` through `nG` and change `Dr`, but that use
should remain within the range supported by laboratory validation. A
different sand requires a new calibration rather than only a change in `Dr`.

The parameter roles above describe their primary influence. Several responses
are coupled, so changing one parameter can affect stiffness, hysteresis,
effective stress, and cyclic strain accumulation simultaneously. Do not treat
these values as independent curve-fitting knobs.

## Fixed reference quantities

The OpenSees material also contains fixed quantities from the Ottawa F65
reference library. They are not arguments of the `RIVASand` command.

| Quantity | Reference value | Purpose |
|---|---:|---|
| Reference Young's modulus | 127339.75550887753 kPa | Sets the reference elastic stiffness |
| Poisson's ratio | 0.30 | Defines the elastic bulk/shear relationship |
| Reference pressure | 101.3 kPa | Pressure used to normalize stiffness |
| Hardening pressure exponent | 0.35 | Controls confinement dependence of hardening |
| Constitutive pressure floor | 0.001 kPa | Prevents evaluation at zero or negative mean effective pressure |
| Default tangent pressure floor | 0.5065 kPa | Regularizes only the tangent returned to OpenSees |

Additional internal cyclic controls remain fixed at their verified Ottawa F65
values. They are deliberately not exposed by the Tcl command. This guide does
not enumerate or formulate those unpublished internal relationships.

## Reference command in kPa units

The following command uses kPa as the stress unit and initializes the material
directly at a material point. OpenSees uses tension-positive stress, so
compressive normal stresses are negative.

```tcl
set matTag 8001
set Dr 0.662962962962963

nDMaterial RIVASand $matTag \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 \
    -nSub 1 -stressScale 1.0 \
    -initialStress -19.4 -19.4 -40.0 0.0 0.0 0.0
```

The `-initialStress` order is

```text
sxx syy szz sxy syz sxz
```

and every stress must use the same unit as the rest of the OpenSees model.
Supplying `-initialStress` without `-stage` selects stage 1 automatically.

## Gravity followed by dynamic loading

For an element model with a gravity/geostatic analysis, create RIVA-Sand in
stage 0 and do not supply `-initialStress`:

```tcl
set matTag 8001
set Dr 0.662962962962963

nDMaterial RIVASand $matTag \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 \
    -rho $rho -nSub 1 -stressScale 1.0 -stage 0
```

Then:

1. Construct the soil column and apply gravity.
2. Complete the gravity analysis successfully so its effective stresses are
   committed at every integration point.
3. Hold the gravity loads constant as required by the analysis model.
4. Activate the nonlinear cyclic stage:

```tcl
updateMaterialStage -material $matTag -stage 1
```

5. Start the dynamic analysis.

Every integration-point copy that uses the material tag initializes from its
own committed effective stress. The committed shear stress becomes the
initial static shear bias. Do not activate stage 1 before a compressive
effective-stress state has been established.

## Optional arguments

| Option | Default | Meaning and use |
|---|---:|---|
| `-rho value` | 0 | Mass density returned by the material; use units consistent with the entire OpenSees model |
| `-nSub integer` | 1 | Fixed constitutive substeps per committed OpenSees increment; start with 1 and increase only after a substep-sensitivity study |
| `-stressScale value` | 1 | Number of analysis stress units per kPa; use 1 for kPa and 1000 for Pa |
| `-pMin value` | 0.001 kPa multiplied by `stressScale` | Constitutive minimum mean effective pressure; changing it changes the material response and requires revalidation |
| `-tangentPMin value` | 0.5065 kPa multiplied by `stressScale` | Minimum pressure used only for the tangent returned to OpenSees; it does not change the constitutive stress update |
| `-stage 0|1` | 0 | Stage 0 is elastic gravity/geostatic setup; stage 1 activates RIVA-Sand cyclic behavior |
| `-initialStress ...` | none | Direct effective-stress initialization in the order `sxx syy szz sxy syz sxz`; if supplied without `-stage`, stage 1 is selected |

The legacy spelling `-noSubsteps` is accepted as an alias for `-nSub`, but new
input files should use `-nSub`.

### Stress-unit examples

- If the analysis uses kPa, use `-stressScale 1.0`. The default `-pMin` is
  0.001 kPa and the default tangent floor is 0.5065 kPa.
- If the analysis uses Pa, use `-stressScale 1000.0`. The corresponding
  defaults are 1 Pa and 506.5 Pa.
- In another stress unit, `-stressScale` is the number of that unit equal to
  1 kPa.

`-stressScale` does not define the OpenSees unit system. OpenSees receives
unitless floating-point values; the user remains responsible for a consistent
force-length-time unit system.

## Constitutive substeps

`-nSub` divides one committed OpenSees strain increment into equal fixed
constitutive subincrements. It does not change the global OpenSees timestep.

Use this procedure:

1. Run the analysis with `-nSub 1`.
2. Repeat representative strong-motion cases with 2 and 4 substeps.
3. Compare peak strain, effective stress, pore pressure, and important station
   responses.
4. Retain the smallest value that gives adequately converged results.

More substeps increase constitutive cost approximately in proportion to
`nSub`; they do not guarantee convergence of a global implicit solution.

## Compatible OpenSees elements

`RIVASand` is a three-dimensional material and can be copied using the
`ThreeDimensional` material type. The implementation is intended for 3D solid
and `u-p` elements that accept a three-dimensional effective-stress material,
including:

- `BrickUP`
- `BBarBrickUP`
- `SSPbrickUP`

In a coupled `u-p` analysis, RIVA-Sand returns the effective skeleton stress.
The element owns the fluid-pressure degree of freedom, permeability, fluid
mass balance, and displacement-pressure coupling.

For station output, compute the pore-pressure ratio from the solved excess
nodal pore pressure and the committed end-of-gravity vertical effective stress.
Do not write the material response `effectivePressureRatio` into the element's
pore-pressure degree of freedom. That material response is only a standalone
skeleton pressure-loss diagnostic.

## Material responses and recorders

The material accepts the following response names:

| Response | Description |
|---|---|
| `stress` | Six effective-stress components in OpenSees order |
| `strain` | Six strain components in OpenSees order; shear entries are engineering shear strains |
| `state` | Complete flattened internal state, primarily for restart and implementation verification |
| `voidRatio` | Current material void ratio |
| `effectivePressureRatio` | Standalone skeleton pressure-loss diagnostic; not the coupled `u-p` pore pressure ratio |
| `reversals` | Accepted constitutive reversal count |
| `compatibilityResidual` | Numerical check on volumetric state compatibility |
| `pressureFloor` | Constitutive pressure floor in the current stress unit |
| `tangentPressureFloor` | Tangent-only pressure floor in the current stress unit |
| `stage` | Current material stage, 0 or 1 |

For example, the response of integration point 1 in an element can be queried
with commands of the following form:

```tcl
set stress [eleResponse $eleTag material 1 stress]
set voidRatio [lindex [eleResponse $eleTag material 1 voidRatio] 0]
set reversals [lindex [eleResponse $eleTag material 1 reversals] 0]
```

The material-point example in `RIVASand_material_point.tcl` demonstrates these
responses through a homogeneous 3D brick.

## Tangent and convergence considerations

During stage 1, RIVA-Sand returns its current pressure- and state-dependent
elastic tangent rather than a consistent algorithmic elastoplastic tangent.
Consequently, strongly nonlinear implicit analyses may require smaller global
steps and may not exhibit quadratic Newton convergence.

At very low confinement, avoid combining a nearly vanishing tangent with an
unnecessarily large penalty constraint. Prefer the `Transformation` constraint
handler when the model topology permits it. If a penalty method is necessary,
scale it relative to the assembled physical stiffness rather than assigning a
very large value automatically.

The two pressure-floor options have different purposes:

- `-tangentPMin` regularizes the tangent supplied to the global solver without
  changing the constitutive stress path.
- `-pMin` prevents the constitutive state from falling below the selected
  residual confinement and therefore changes the predicted response.

A site-specific change to `-pMin` is a model variant and should be checked
against material-point and site-response benchmarks.

## Common input mistakes

- Writing `Dr=66.3` instead of `Dr=0.663`.
- Using `-stressScale 1` while all stresses are entered in Pa.
- Entering compressive initial normal stresses as positive values.
- Selecting `-stage 1` without supplying `-initialStress`.
- Starting dynamics without calling `updateMaterialStage` after gravity.
- Assuming the other ten required parameters change automatically with `Dr`.
- Treating `effectivePressureRatio` as the nodal pore pressure in a `u-p`
  element.
- Increasing `-pMin` for numerical convenience without rechecking the physical
  response.
- Assuming the Ottawa F65 row is calibrated for another sand.

## Calibration limitations

The parameter row in this guide is calibrated for Ottawa F65 sand using the
available cyclic direct simple shear database. It is not a universal sand
parameter set. The strongest existing evidence concerns medium-dense to dense
conditions covered by the calibration and validation tests. Known limitations
include imperfect transfer of cyclic resistance across some dense,
statically-biased cases and underprediction of irreversible pore-pressure
buildup in some biased histories.

For a different material, calibration should use, at minimum, relative-density
and void-ratio information, small-strain stiffness or shear-wave velocity,
modulus-reduction data, and cyclic laboratory histories covering the intended
confinement, density, cyclic amplitude, and static-bias range.

## Verification example

From the OpenSees source root, run the supplied material-point example with a
build that contains RIVA-Sand:

```bash
./build-rivasand/OpenSees EXAMPLES/RIVASand/RIVASand_material_point.tcl
```

The example writes `RIVASand_material_point.csv`. The verification files in
`EXAMPLES/RIVASand/golden_data/` can then be used to confirm that the native
OpenSees response matches the frozen reference histories.
