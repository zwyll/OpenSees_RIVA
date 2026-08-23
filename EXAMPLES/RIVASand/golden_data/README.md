# Golden verification data

These files are deterministic material-point results from the standalone RIVA-Sand
kernel with its host-step reversal detector enabled.  They are intended to
verify a native Hercules implementation increment by increment.

## Files

- `manifest.json`: case definitions, row counts, summaries, and SHA-256 hashes.
- `reference_parameters.json`: all 112 frozen constitutive parameters and four
  port-only event controls.
- `state_schema.json`: all 38 state fields, types, tensor component order, and
  the symmetric-tensor double-contraction rule.
- `cyclic_zero_bias_reference.csv`: zero-bias cyclic baseline.
- `cyclic_bias025_reference.csv`: reference-density biased response.
- `cyclic_bias0375_dense.csv`: dense biased response.
- `cyclic_bias025_four_substeps.csv`: fixed-substep and reversal-event case.
- `mixed_3d_bias015.csv`: nonproportional three-dimensional tensor path.
- `RIVASAND_VERIFICATION_INDEX.xlsx`: human-readable case, parameter, and
  state inventory; the CSV/JSON files remain the machine-readable authority.

Every CSV starts with row `step=0`, the activated initial state.  Subsequent
rows contain one Hercules/global strain increment and the expected state after
that increment.  Tensor components are true tensor values in this order:

```text
xx, yy, zz, xy, yz, xz
```

The `deps_*` fields therefore contain tensor shear strain, not engineering
shear.  For example, `deps_xz = Delta gamma_xz/2`.

## Acceptance sequence

1. First reproduce `step=0` exactly enough to confirm signs, K0 stress, static
   bias, and initialization.
2. Compare stress and event counters with `--mode core`.
3. Compare all history variables with `--mode full`.
4. Use tight binary64 tolerances while porting.  Relax tolerances only after a
   documented comparison of operation ordering and compiler math behavior.

Example:

```bash
python3 verify_hercules_output.py golden_data hercules_results --mode core
python3 verify_hercules_output.py golden_data hercules_results --mode full
```

Regenerate only from the frozen package root:

```bash
python3 -m rivasand_port.generate_golden_data
```

If regenerated hashes change, treat that as a model-source change and review it
before updating any Hercules baseline.
