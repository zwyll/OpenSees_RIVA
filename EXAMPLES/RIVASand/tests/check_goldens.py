"""Read-only golden replay check: run GoldenTrace on each golden CSV and
report the max relative drift of every kernel-derived column vs the stored
values. PASS if drift <= 1e-9 for every case (bit-preservation gate).

Usage: python check_goldens.py <golden_dir> <GoldenTrace.exe>
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

GD = Path(sys.argv[1])
TRACE = sys.argv[2]

schema = json.loads((GD / "state_schema.json").read_text())
flat = []
for f in schema["fields"]:
    if f["components"]:
        flat += [f"{f['name']}_{c}" for c in f["components"]]
    else:
        flat.append(f["name"])

INPUT_COLS = {"case_id", "step", "time_s", "path_coordinate", "nsub",
              "engineering_gamma_xz"}
manifest = json.loads((GD / "manifest.json").read_text())
worst_all = 0.0
fails = 0
for case in manifest["cases"]:
    src = GD / case["file"]
    with open(src, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        old_rows = [dict(zip(header, r)) for r in rd]
    out = subprocess.run([TRACE, str(src)], capture_output=True, text=True)
    if out.returncode != 0:
        print(f"TRACE FAILED {case['file']}: {out.stderr[:300]}")
        fails += 1
        continue
    trows = [ln.split(",") for ln in out.stdout.strip().splitlines()]
    assert len(trows) == len(old_rows), (case["file"], len(trows),
                                         len(old_rows))
    kcols = [c for c in header if c not in INPUT_COLS]
    worst = 0.0
    wcol = ""
    for old, tr in zip(old_rows, trows):
        # trace output: info triple + flat state values, same order as regen
        tvals = dict(zip(["accepted_substeps", "reversal_registered",
                          "compatibility_residual"] + flat, tr))
        for c in kcols:
            if c not in tvals:
                continue
            a, b = float(old[c]), float(tvals[c])
            d = abs(a - b) / (1.0 + abs(a))
            if d > worst:
                worst, wcol = d, c
    worst_all = max(worst_all, worst)
    ok = worst <= 1e-9
    fails += 0 if ok else 1
    print(f"{'PASS' if ok else 'FAIL'} {case['file']}: max drift "
          f"{worst:.3e} ({wcol})")
print(f"{'ALL PASS' if fails == 0 else str(fails) + ' FAILURES'} "
      f"(worst {worst_all:.3e})")
sys.exit(1 if fails else 0)
