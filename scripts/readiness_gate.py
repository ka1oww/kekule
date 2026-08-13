"""Run every local check required for the narrow Kekule preview."""

import os
import subprocess
import sys


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COMMANDS = (
    ("compile", [sys.executable, "-m", "compileall", "-q", "."]),
    ("66-case invariants", [sys.executable, "tests/invariants.py"]),
    (
        "32 validation tests",
        [
            sys.executable,
            "-m",
            "unittest",
            "discover",
            "-s",
            "tests",
            "-p",
            "test_validation.py",
            "-v",
        ],
    ),
    ("preview manifest", [sys.executable, "tests/run_preview_manifest.py"]),
)


def main():
    for name, command in COMMANDS:
        print(f"\n== {name} ==", flush=True)
        subprocess.run(command, cwd=ROOT, check=True)
    print("\nKekule readiness gate passed.")


if __name__ == "__main__":
    main()
