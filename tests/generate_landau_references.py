#!/usr/bin/env python3
"""
Script to automatically generate reference files for regression tests for the Landau free energy.
"""

import shutil
import sys
from pathlib import Path
import tempfile

from test_landau import SOLVERS, INTEGRATORS, ROOT, REF, run_simulation

def main():
    # Ensure binary directory exists
    if not (ROOT / "build" / "bin").exists():
        print("Error: Build directory not found. Please run CMake configure and build first.")
        print(f"  Expected: {ROOT / 'build' / 'bin'}")
        sys.exit(1)

    # Ensure all binaries exist
    missing_bins = [name for name, path in SOLVERS.items() if not path.exists()]
    if missing_bins:
        print(f"Error: Missing binary(ies): {', '.join(missing_bins)}")
        print(f"  Please build the project first: cd {ROOT / 'build'} && cmake --build .")
        sys.exit(1)

    # Ensure reference directory exists
    REF.mkdir(parents=True, exist_ok=True)

    print(f"Generating reference files for regression tests for the Landau free energy...")
    print(f"Reference directory: {REF}")
    print()

    total = len(SOLVERS) * len(INTEGRATORS)
    success = 0

    # Run all simulation combinations
    for solver_key, integrator in [(k, i) for k in SOLVERS.keys() for i in INTEGRATORS]:
        ref_file = REF / f"landau_{integrator}_{solver_key}.dat"
        solver = SOLVERS[solver_key]

        print(f"Generating {ref_file.name}...", end=" ", flush=True)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            result = run_simulation(tmp_path, solver, integrator)

            if result and result.exists():
                # Copy to reference directory
                shutil.copy(result, ref_file)
                print("Simulation completed, reference file generated.")
                success += 1
            else:
                print("Error: Simulation did not produce expected output.")

    print()
    print(f"Generated {success}/{total} reference files")

    if success == total:
        print("All reference files generated successfully!")
        return 0
    else:
        print(f"Warning: {total - success} reference file(s) failed to generate")
        return 1


if __name__ == "__main__":
    sys.exit(main())
