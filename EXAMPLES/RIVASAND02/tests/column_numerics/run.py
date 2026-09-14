"""Run the SSPbrickUP and branch-reversal research regression checks.

Requires Python with NumPy and a built OpenSees executable. Supply --probe,
or --cmake-build for a Unix Makefiles build using the same compiler as OpenSees.
The probe only exposes existing Domain/Element operations to the test scripts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def build_probe(build, output):
    flags = {}
    for line in (build / "CMakeFiles/OPS_InterpTcl.dir/flags.make").read_text().splitlines():
        if " = " in line:
            key, value = line.split(" = ", 1)
            flags[key] = shlex.split(value)
    link = shlex.split((build / "CMakeFiles/OpenSees.dir/link.txt").read_text())
    command = [link[0], "-std=c++17", "-O2", "-fPIC"]
    command += flags["CXX_DEFINES"] + flags["CXX_INCLUDES"] + flags["CXX_FLAGS"]
    if sys.platform == "darwin":
        command += ["-dynamiclib", "-undefined", "dynamic_lookup"]
    elif sys.platform.startswith("linux"):
        command += ["-shared"]
    else:
        raise RuntimeError("Build the Tcl probe separately and supply --probe on this platform")
    command += [flag for flag in ("-static-libgcc", "-static-libstdc++") if flag in link]
    library = output / "columnprobe.so"
    command += [str(HERE / "column_probe.cpp"), "-o", str(library)]
    subprocess.run(command, cwd=build, check=True)
    return library


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opensees", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--probe", type=Path)
    source.add_argument("--cmake-build", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    executable = args.opensees.resolve()
    probe = args.probe.resolve() if args.probe else build_probe(args.cmake_build.resolve(), output)
    checks = []
    for script, destination, marker, extra in [
        ("probe_column_contracts.tcl", "contracts", "CONTRACT_PROBES_COMPLETE", ["strict"]),
        ("rotating_path.tcl", "rotating_results", "ROTATING_PATH_CHECKS_COMPLETE", []),
        ("rotation_covariance.tcl", "rotation_results", "ROTATION_COVARIANCE_COMPLETE", []),
        ("guard_checks.tcl", "guard_results", "PASS: research latch rejection", []),
    ]:
        result = subprocess.run(
            [str(executable), str(HERE / script), str(probe), str(output / destination), *extra],
            cwd=output, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180,
        )
        (output / (Path(script).stem + ".log")).write_text(result.stdout)
        if result.returncode != 0 or marker not in result.stdout:
            raise RuntimeError(f"{script} failed:\n{result.stdout[-5000:]}")
        checks.append(script)
        print(f"PASS: {script}", flush=True)
    subprocess.run([sys.executable, str(HERE / "verify.py"), str(output)], check=True)
    report = {
        "status": "passed",
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "checks": checks,
        "scope": "Element contracts, prescribed paths and research adapter guards; no production acceptance",
    }
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
