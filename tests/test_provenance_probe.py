#!/usr/bin/env python3
"""Offline regression tests using synthetic fixtures only."""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROBE = os.path.join(ROOT, "provenance_probe.py")
EXAMPLES = os.path.join(ROOT, "examples")
fails = []


def run(args):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        output = tmp.name
    result = subprocess.run(
        [sys.executable, PROBE] + args + ["--json", output],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"exit {result.returncode}\n{result.stdout}\n{result.stderr}"
        )
    with open(output) as handle:
        report = json.load(handle)
    os.unlink(output)
    return report, result.stdout


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        fails.append(name)


target_original = os.path.join(EXAMPLES, "target_original.jsonl")
target_substituted = os.path.join(EXAMPLES, "target_substituted.jsonl")
clean_original = os.path.join(EXAMPLES, "clean_original.jsonl")
clean_substituted = os.path.join(EXAMPLES, "clean_substituted.jsonl")
screen_demo = os.path.join(EXAMPLES, "screen_demo.jsonl")

print("test: synthetic target-dependent construction")
report, output = run(
    ["intervene", "--original", target_original, "--substituted", target_substituted]
)
check("all channel values move", report["channel_moved_frac"] == 1.0)
check("geometry is checked", report["geometry_checked"])
check("target dependence is detected", report["verdict"] == "TARGET_DEPENDENCE_DETECTED")
check("comparison is exact", report["exact"] is True)

print("\ntest: synthetic design-valid construction")
report, output = run(
    ["intervene", "--original", clean_original, "--substituted", clean_substituted]
)
check("clean channel is unchanged", report["channel_moved_frac"] == 0.0)
check("unchanged is not over-certified", report["verdict"] == "NO_DEPENDENCE_DETECTED")
check("invariance caveat is printed", "No finite unchanged substitution certifies" in output)

print("\ntest: observational screen")
report, output = run(["screen", "--train", screen_demo, "--test", screen_demo])
check("screen declares itself observational", "OBSERVATIONAL" in output)
check("screen issues no provenance verdict", "NO PROVENANCE VERDICT IS ISSUED" in output)
check("screen does not print target-derived", "TARGET-DERIVED" not in output)
check("screen does not print design-valid", "DESIGN-VALID" not in output)
check("screen report is inexact", report["exact"] is False)

print("\ntest: documented CLI")
for command in (
    ["intervene", "--original", target_original, "--substituted", target_substituted],
    ["screen", "--train", screen_demo, "--test", screen_demo],
):
    result = subprocess.run([sys.executable, PROBE] + command, capture_output=True)
    check(f"{command[0]} runs", result.returncode == 0)
result = subprocess.run([sys.executable, PROBE, "--version"], capture_output=True, text=True)
check("version option works", result.returncode == 0 and "provenance_probe" in result.stdout)

print("\ntest: malformed inputs and adversarial controls")
with tempfile.TemporaryDirectory() as directory:
    left = os.path.join(directory, "left.jsonl")
    right = os.path.join(directory, "right.jsonl")

    def write(path, rows):
        with open(path, "w") as handle:
            for row in rows:
                handle.write(json.dumps(row) + "\n")

    def compare(left_rows, right_rows):
        write(left, left_rows)
        write(right, right_rows)
        return run(["intervene", "--original", left, "--substituted", right])[0]

    base = {"chain": "rare", "channel": [0] * 200, "geom": [[0, 0, 0]] * 200}
    rare = dict(base, channel=[1] + [0] * 199)
    report = compare([base], [rare])
    check("0.5 percent movement is detected", report["n_changed_values"] == 1)

    invariant = {"chain": "class", "channel": [0, 0, 1, 1], "geom": [0] * 4}
    report = compare([invariant], [invariant])
    check("class-preserving change cannot certify validity", report["verdict"] == "NO_DEPENDENCE_DETECTED")

    missing = {"chain": "missing", "channel": [0]}
    report = compare([base, missing], [rare, dict(missing, channel=[1])])
    check("partial geometry is inconclusive", report["verdict"] == "INCONCLUSIVE_UNVERIFIED_CONTROL")
    report = compare([missing], [dict(missing, channel=[1])])
    check("missing geometry is inconclusive", report["verdict"] == "INCONCLUSIVE_UNVERIFIED_CONTROL")
    changed_geometry = dict(rare, geom=[[1, 0, 0]] + [[0, 0, 0]] * 199)
    report = compare([base], [changed_geometry])
    check("geometry movement invalidates intervention", report["verdict"] == "INVALID_INTERVENTION")

    invalid_cases = [
        ("empty input", [], []),
        ("duplicate IDs", [base, base], [rare, rare]),
        ("empty channels", [dict(base, channel=[])], [dict(base, channel=[])]),
        ("chain ID mismatch", [base], [dict(rare, chain="other")]),
        ("channel length mismatch", [base], [dict(rare, channel=[1])]),
        ("geometry length mismatch", [base], [dict(rare, geom=[[0, 0, 0]])]),
    ]
    for name, left_rows, right_rows in invalid_cases:
        write(left, left_rows)
        write(right, right_rows)
        result = subprocess.run(
            [sys.executable, PROBE, "intervene", "--original", left, "--substituted", right],
            capture_output=True,
        )
        check(f"{name} is rejected", result.returncode != 0)

print(f"\n{'ALL 27 TESTS PASSED' if not fails else str(len(fails)) + ' FAILURE(S)'}")
sys.exit(1 if fails else 0)
