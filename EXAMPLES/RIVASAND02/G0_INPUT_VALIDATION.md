# G0 input migration and verification — 2026-09-08

Both `RIVASand` and `RIVASAND02` now require this sequence after the tag:

```text
Dr G0 M kd h m zeta eMax eMin Q R nG
```

`G0` is positive and dimensionless. The reference value
`483.48301127222084` retains the previous stiffness exactly, corresponding to
a reference shear modulus of `48976.82904187597 kPa` at `101.3 kPa` mean
effective pressure. The same G0 is used in Pa and kPa analyses. Changing it
scales elastic shear and bulk stiffness together without changing the pressure
exponent, pressure guards, cyclic controls, or integration event rules.

## Compatibility

Insert G0 after Dr in every previous eleven-value command. Existing optional
arguments retain their meanings. The original command is case-sensitive
`RIVASand`; the successor is `RIVASAND02`.

New serialization contains 135 values for the original adapter and 183 for
the successor, including G0 and an adapter-revision marker. Pre-G0 databases
are incompatible: rerun initialization. Constitutive state sizes are unchanged
(93 and 139 flattened values respectively). The successor's kernel revision
remains 4; the new input/serialization marker is separate.

The successor now also derives its reference-density cache from the supplied
eMax/eMin, matching the existing Hercules correction and private reference.
Reference density bounds preserve the accepted histories. Results obtained
with non-reference bounds before this fix must be rerun.

## Completed checks

Built and tested with GNU C++ 15 on macOS:

- Full OpenSees executable rebuild, including both material adapters.
- All five original and six successor native golden histories, without
  modifying reference data or changing their existing tolerances.
- Native state, geostatic-admission, plastic-activity-gate, and option tests.
- Original and successor stage activation, plus successor database restart,
  no-bias-volume, and all three BiasVolume-mode Tcl regressions.
- The new `tests/RIVASand_G0.tcl` tests both generations at G0=350, reference,
  and 700 in Pa and kPa: analytical elastic response, material isolation,
  element copies, nonlinear save/restore, and rejection of invalid G0 and
  obsolete commands. It additionally checks four eMax/eMin combinations and
  retention of the material-specific reference density after restart.

Run the input regression with the newly built executable from a scratch
working directory (the test creates database and response files):

```sh
/path/to/OpenSees /path/to/EXAMPLES/RIVASAND02/tests/RIVASand_G0.tcl
```

The added input uses a setup-time stiffness cache; no new arithmetic is added
inside the constitutive substep. These checks establish input and material-point
behavior, not full LEAP or 3-D boundary-value qualification. Non-reference G0
values still require physical calibration and renewed field-response checks;
this change does not demonstrate reduced accelerations or unchanged pore pressure.

See [the user guide](../RIVASand/RIVASand_USER_GUIDE.md) for parameter meanings,
units, and calibrated values.
