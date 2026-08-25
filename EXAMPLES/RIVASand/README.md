# RIVA-Sand for OpenSees

**RIVA-Sand** denotes the **R**eversible--**I**rreversible **V**olumetric
**A**ccumulation model for Sand. `RIVASand` is its three-dimensional,
effective-stress `NDMaterial` adapter. The constitutive update is the
same allocation-free, fixed-substep Forward Euler kernel used by the verified
Hercules port. OpenSees-specific code supplies engineering-shear conversion,
trial/committed state management, gravity-to-dynamic activation, recorders,
and channel serialization.

## Command

```tcl
nDMaterial RIVASand tag Dr M kd h m zeta eMax eMin Q R nG \
    <-rho value> <-nSub integer> <-stressScale value> \
    <-pMin value> <-tangentPMin value> <-stage 0|1> \
    <-initialStress sxx syy szz sxy syz sxz>
```

The reference Ottawa F65 row is

```text
Dr M kd h m zeta eMax eMin Q R nG
0.662962962962963 1.25 1.125 122.44207260469 0.945 0.025 0.78 0.51 10 1.5 0.65
```

`Dr` is a fraction in `[0,1]`, not percent. The initial internal void ratio is

```text
e0 = eMax - Dr (eMax-eMin).
```

The normal OpenSees sign convention is used: tensile normal stress is positive
and a compressive `-initialStress` is negative. The strain vector is
`[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]`; the adapter divides
the engineering shear components by two before calling the tensor kernel.

Options:

- `-rho`: mass density returned to the element. Default `0`; use the density
  convention required by the selected OpenSees element.
- `-nSub`: fixed constitutive substeps per committed OpenSees increment.
  Default `1`. The reversal detector runs once per host increment, not once per
  substep.
- `-stressScale`: converts the frozen kPa calibration to the model's stress
  unit. Use `1` for kPa and `1000` for Pa. Default `1`.
- `-pMin`: optional **constitutive** minimum mean effective pressure, in the
  current stress unit. It is used by the explicit stress update and therefore
  changes the response. The frozen V8 default remains `0.001 kPa` times
  `stressScale`.
- `-tangentPMin`: minimum pressure used only to form the elastic tangent
  returned to OpenSees. It does not change the V8 stress update. The default is
  `p_ref/200` (`0.5065 kPa` when `stressScale=1`) and is never allowed below
  `-pMin`.
- `-stage`: `0` for elastic gravity/geostatic setup and `1` for V8 dynamics.
  Default `0`.
- `-initialStress`: direct material-point initialization in OpenSees stress
  order. Supplying it without `-stage` selects stage 1.

`-stage 1` at construction requires `-initialStress`. In an element model,
start at stage 0, finish and commit gravity, and then use

```tcl
updateMaterialStage -material $matTag -stage 1
```

Each integration-point copy then initializes from its own committed effective
stress. Its committed normal stress is retained as the geostatic reference;
committed shear stress is interpreted as the static bias.

The OpenSees `MaterialStageParameter` implementation must scan all elements,
not stop after the first match. `tests/RIVASand_stage_activation.tcl`
uses two elements sharing one material tag and fails if either copy remains at
stage 0.

## u-p elements

The material returns **effective skeleton stress** and a 3D tangent, so it can
be copied by OpenSees `BrickUP`, `BBarBrickUP`, and `SSPbrickUP` elements through
`getCopy("ThreeDimensional")`. The element continues to own its pressure DOF,
fluid mass balance, permeability, and displacement-pressure coupling.

Do not write the material response `effectivePressureRatio` into the element's
pore-pressure DOF. It is only the standalone skeleton diagnostic
`1-p'/p'_anchor`. For a coupled u-p station output, calculate
`r_u = Delta u/sigma'_v0` from the solved nodal pore pressure and the committed
end-of-gravity vertical effective-stress reference.

## Tangent and analysis type

V8 is integrated explicitly inside the material. Stage 1 returns its current
pressure/state-dependent elastic tangent, not a consistent algorithmic
elastoplastic tangent. This is suitable for explicit or transient workflows
where the constitutive substep is the intended integration scheme. A strongly
nonlinear implicit static solve can require smaller global steps and may not
show quadratic Newton convergence.

For low-confinement u-p analyses, do not combine a nearly vanishing tangent
with an unnecessarily large penalty constraint. Prefer the Transformation
constraint handler where the model topology permits it. If Penalty is needed,
choose its scale relative to the assembled physical stiffness rather than
using `1e18` automatically. A site-specific `-tangentPMin 2.0` can regularize
the OpenSees tangent without silently replacing the frozen constitutive
pressure floor. If the explicit update itself must retain residual confinement,
use `-pMin 2.0` as a clearly identified site-response model variant and
revalidate its low-pressure/post-liquefaction response; a tangent-only floor
cannot prevent the constitutive state from reaching the frozen 0.001 kPa apex.

## Responses

The material accepts these response names:

- `stress`, `strain`
- `state` (the 93 flattened state values in formulation order)
- `voidRatio`
- `effectivePressureRatio` (standalone skeleton diagnostic only)
- `reversals`
- `compatibilityResidual`
- `pressureFloor` (constitutive `p_min`)
- `tangentPressureFloor` (OpenSees tangent-only floor)
- `stage`

## Files and verification

- `RIVASand_material_point.tcl` drives the zero-bias, two-cycle reference
  history through a homogeneous `bbarBrick`, thereby testing the production
  OpenSees element/material interface as well as the adapter.
- `RIVASand_FORMULATION.md` contains every equation, parameter, state
  variable, update order, and initialization rule.
- `tests/RIVASandKernelStateTest.cpp` checks reference/custom equivalence,
  active `zeta`, trial-state isolation, restart, and zero increments.
- `tests/RIVASandGoldenReplay.cpp` replays all five frozen 1D/3D golden paths,
  including biased, dense, four-substep, and nonproportional histories.
- `tests/RIVASand_stage_activation.tcl` verifies that two `SSPbrickUP`
  elements accept the material, and that `updateMaterialStage` reaches both
  element-held copies and initializes each from committed stage-0 stress.
- `golden_data/` is a self-contained copy of the machine-readable frozen
  oracle, parameter inventory, and state schema.

The dependency-free kernel tests can be built from the OpenSees source root:

```bash
c++ -std=c++17 -O2 EXAMPLES/RIVASand/tests/RIVASandKernelStateTest.cpp \
  -o RIVASandKernelStateTest
c++ -std=c++17 -O2 EXAMPLES/RIVASand/tests/RIVASandGoldenReplay.cpp \
  -o RIVASandGoldenReplay
./RIVASandKernelStateTest
./RIVASandGoldenReplay
```

The Tcl driver has also been run with a fully linked OpenSees 3.8.0 executable.
Its 65 output states match `cyclic_zero_bias_reference.csv` with maximum
absolute differences of `3.2e-8` kPa in shear stress, `2.9e-6` kPa in normal
stress, and `1.1e-7` in the skeleton pressure-loss diagnostic. Reversal counts
match exactly and the compatibility residual differs by less than `2e-19`.

The kernel calibration SHA-256 is
`9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83`.
