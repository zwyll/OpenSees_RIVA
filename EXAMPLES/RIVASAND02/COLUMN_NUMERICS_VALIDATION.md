# Column numerical corrections

The coupled free-field investigation identified two SSPbrickUP implementation
errors and a reversal-detection sensitivity in RIVASAND02. This branch fixes
the element errors and adds the separately named built-in research material
`RIVASAND02BranchReversalResearch`.

## SSPbrickUP corrections

- `update()` returns the material's `setTrialStrain()` result and reports the
  element tag when the material rejects a trial. Previously it always returned
  success. This detects a rejected numerical update; ordinary yielding or
  liquefaction is not itself an error condition.
- The `betaKcomm` Rayleigh contribution uses the solid-displacement block of
  the last committed element stiffness. Previously it used current trial
  stiffness. The displacement indices skip the fourth, fluid degree of freedom
  at each node. This correction is limited to `betaKcomm`; the other Rayleigh
  contributions and the hydraulic/coupling blocks retain their existing behavior.

The focused regression makes an uncommitted strain change after saving the
committed stiffness and compares the requested damping block against that
saved matrix. It also injects an invalid strain and requires Domain::update
to report failure. The original executable gives a damping relative error
of 1.5183 and incorrectly accepts the invalid trial. Both checks pass with
the corrections (zero damping error, invalid trial rejected).

## Research material usage

Replace the material identifier with `RIVASAND02BranchReversalResearch` to
select the corrected event reference. It is built into the same OpenSees
material source; no external material plugin is needed. It has a distinct
class tag, and material copies retain its identity.

For example, this creates the dense-layer research material at stage 0,
using kPa stress units and the fixed inputs from the numerical study:

```tcl
nDMaterial RIVASAND02BranchReversalResearch 101 \
    0.90 900 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 4 4 0.65 \
    -rho 2.073519843851659 -nSub 10 -TanType 0 -noBiasVolume
```

Follow the existing gravity and material-stage procedure. The twelve physical
parameters, units, recorders, and tangent selector have the same meanings as
in RIVASAND02. The corrected reference follows accumulated deviatoric strain
along a loading branch, so gradually rotating increment directions can still
register reversal. Physical parameters are unchanged.

This changes constitutive event handling and remains a **research successor**.
The original `RIVASAND02` command retains its previous event rule and restart
format. The research command rejects `-reversalLatch`. Its `sendSelf` and
`recvSelf` reject database save/restore and channel transfer, including
partitioned execution that needs material serialization: six historical state
slots now have different meanings. Start research cases from initialization.

## Verification scope

Two prescribed shear cycles with simultaneous vertical strain register four
reversals at 128, 512, 2,048, and 8,192 steps per cycle with the correction.
The previous material registers four at 128 and zero at the finer resolutions.
Pure-shear stress, pressure, and reversal histories are identical between the
two materials at each resolution. Rotating the combined path for dense and
loose rows changes the transformed stress by less than 1.5e-13 kPa in the
tested cases. Speculative uncommitted trials do not change accepted results;
copy identity, latch rejection, and checkpoint rejection are also checked.

The earlier fixed-input, 20-element, 70-second column study supports
`NewmarkExplicit 0.5`, `Linear`, `Transient`, **dt = 0.0025 s**, **nSub = 10**,
`-TanType 0`, and `-noBiasVolume` for working runs in that specific setup.
Relative to the 0.0003125 s benchmark, maximum PGA and pressure-peak changes
are 0.0262% and 0.0615%; maximum spectral log RMS is 0.001139.
The 0.005 s explicit case fails with nSub 10, 20, 40, 80, 160, and 320.
These time steps are case-specific and must be checked again after stiffness,
mesh, or other relevant inputs change.

The built-in material was also replayed through the complete 70-second column
at 0.0025 s and nSub = 10. All six recorder histories (absolute acceleration,
surface acceleration, displacement, pore pressure, strain, and stress) are
exactly identical to the previously tested plugin at all 28,001 recorded
times. Seven existing Tcl adapter regressions and both standalone kernel
state-contract tests pass. The compact [verification record](COLUMN_NUMERICS_RESULTS.json)
identifies the tested executable and the reference plugin.

The local executable was linked from the existing CMake build after recompiling
the affected element, material, Tcl/generic factories, and active class broker.
The optional Xara runtime broker is not part of this executable; a separate
syntax check cannot run in this checkout because its existing `Hash.h`
dependency is absent. No Xara or distributed-execution acceptance is claimed.

This is numerical verification, not experimental recalibration or production
acceptance. General material validation, mesh convergence, full-duration
implicit convergence, checkpoint support, other bias-volume configurations,
and Hercules GPU/MPI validation remain outstanding. No new h/m fit is included.

## Reproduce the focused checks

Build OpenSees normally with this source. From the repository root, with
Python and NumPy installed:

```sh
python EXAMPLES/RIVASAND02/tests/column_numerics/run.py \
    --opensees build/OpenSees --cmake-build build \
    --output build/column_numerics_checks
```

The probe build helper supports existing Unix Makefiles CMake builds on macOS
and Linux and uses that build's compiler and include flags. Alternatively,
build `column_probe.cpp` as a Tcl extension against the same OpenSees build
and supply `--probe /path/to/columnprobe.so`. This extension only exposes
Domain/Element operations for verification. It does not implement a material.

The runner checks completion markers as well as process status, verifies
finite recorded responses, and writes compact JSON results. It includes the
SSPbrickUP contracts, the prescribed-path refinement and rotation tests, and
research adapter guards. Existing RIVASand/RIVASAND02 state, stage, G0,
restart, bias-volume, and continuum-adapter checks remain applicable.
