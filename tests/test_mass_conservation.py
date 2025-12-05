import numpy as np
import subprocess
import shutil
from pathlib import Path
import itertools
import pytest

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
REF  = DATA / "reference"

BIN_1D = ROOT / "build" / "bin" / "ch_1D"
BIN_2D = ROOT / "build" / "bin" / "ch_2D"

SOLVERS = {
    "1D": BIN_1D,
    "2D": BIN_2D,
}
INTEGRATORS = ["euler", "pseudospectral"]

def load_field(path):
    """Load a binary field written as doubles."""
    return np.fromfile(path, dtype=np.float64)

def run_simulation(tmp_path, solver, integrator):
    """Run simulation with given solver + integrator inside tmp_path."""
    # prepare per-test config
    config = tmp_path / "simple_landau.toml"
    shutil.copy(DATA / "simple_landau.toml", config)

    # modify integrator inside tmp config
    text = config.read_text()
    text = f'\nintegrator = "{integrator}"\n' + text
    config.write_text(text)

    subprocess.run([str(solver), str(config)], cwd=tmp_path, check=True)
    return tmp_path / "last_0.dat"

@pytest.mark.parametrize("solver_key", ["1D", "2D"])
@pytest.mark.parametrize("integrator", INTEGRATORS)
def test_mass_conservation(tmp_path, solver_key, integrator):
    solver = SOLVERS[solver_key]
    final_conf = run_simulation(tmp_path, solver, integrator)

    phi0 = load_field(tmp_path / "init_0.dat")
    phi1 = load_field(tmp_path / "last_0.dat")

    assert abs(phi0.sum() - phi1.sum()) < 1e-12

@pytest.mark.parametrize("solver_key", ["1D", "2D"])
@pytest.mark.parametrize("integrator", INTEGRATORS)
def test_regression_last_configuration(tmp_path, solver_key, integrator):
    solver = SOLVERS[solver_key]
    final_path = run_simulation(tmp_path, solver, integrator)
    new = load_field(final_path)

    ref_file = REF / f"landau_{integrator}_{solver_key}.dat"
    ref = load_field(ref_file)

    np.testing.assert_allclose(new, ref, rtol=1e-12, atol=1e-14)
