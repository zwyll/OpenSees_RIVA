"""Regenerate the golden oracle CSVs from the CURRENT kernel (V9 re-freeze).

Keeps the input schedule (deps/nsub/time) and header of each existing case,
replaces every kernel-derived column, backs up the old files, quantifies the
drift, and updates manifest.json metadata (+sha256 where present).

Usage: python regen_goldens.py <golden_dir> <GoldenTrace.exe>
"""
import csv
import hashlib
import json
import shutil
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
assert len(flat) == 94, len(flat)

INPUT_COLS = {"case_id", "step", "time_s", "path_coordinate", "nsub",
              "engineering_gamma_xz"}
TRACE_INFO = {"accepted_substeps", "reversal_registered",
              "compatibility_residual"}

backup = GD / "v8_backup"
backup.mkdir(exist_ok=True)
manifest = json.loads((GD / "manifest.json").read_text())
report = []

for case in manifest["cases"]:
    fn = case["file"]
    src = GD / fn
    shutil.copy2(src, backup / fn)

    with open(src, newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        old_rows = [dict(zip(header, r)) for r in rd]

    out = subprocess.run([TRACE, str(src)], capture_output=True, text=True)
    if out.returncode != 0:
        print(f"TRACE FAILED for {fn}: {out.stderr}")
        sys.exit(1)
    trows = [ln.split(",") for ln in out.stdout.strip().splitlines()]
    assert len(trows) == len(old_rows), (fn, len(trows), len(old_rows))

    sigv0 = -float(old_rows[0]["stress_zz"])
    p0 = float(old_rows[0]["mean_effective_pressure"])

    drift = 0.0
    ru_max = -9e9
    new_rows = []
    for old, tr in zip(old_rows, trows):
        st = dict(zip(flat, map(float, tr[3:])))
        info = dict(zip(["accepted_substeps", "reversal_registered",
                         "compatibility_residual"], tr[:3]))
        p = -(st["stress_xx"] + st["stress_yy"] + st["stress_zz"]) / 3.0
        sxx, syy, szz = (st["stress_xx"] + p, st["stress_yy"] + p,
                         st["stress_zz"] + p)
        q = (1.5 * (sxx**2 + syy**2 + szz**2
                    + 2 * (st["stress_xy"]**2 + st["stress_yz"]**2
                           + st["stress_xz"]**2))) ** 0.5
        ru_v = 1.0 + st["stress_zz"] / sigv0
        ru_max = max(ru_max, ru_v)
        derived = {"mean_effective_pressure": p, "q": q,
                   "ru_vertical": ru_v, "ru_mean": 1.0 - p / p0}
        row = {}
        for col in header:
            if col in INPUT_COLS or col.startswith(("deps_", "strain_total_")):
                row[col] = old[col]
            elif col in TRACE_INFO:
                row[col] = info[col]
            elif col in derived:
                row[col] = repr(derived[col])
            elif col in st:
                row[col] = repr(st[col])
            else:
                print(f"UNMAPPED column {col} in {fn}")
                sys.exit(1)
            if col not in INPUT_COLS and not col.startswith(
                    ("deps_", "strain_total_")) and col != "case_id":
                try:
                    o, n = float(old[col]), float(row[col])
                    d = abs(n - o) / max(1.0, abs(o))
                    drift = max(drift, d)
                except ValueError:
                    pass
        new_rows.append(row)

    with open(src, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(new_rows)

    case["maximum_ru_vertical"] = ru_max
    case["final_amplitude_reversals"] = int(
        float(new_rows[-1]["amplitude_reversals"]))
    if "sha256" in case:
        case["sha256"] = hashlib.sha256(src.read_bytes()).hexdigest()
    report.append((fn, drift, ru_max))

(GD / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True))
print("V9 re-freeze complete (old files in v8_backup/):")
for fn, drift, ru in report:
    print(f"  {fn:44s} max drift = {drift:.3e}   new ru_max = {ru:.4f}")
