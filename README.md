# Cahn-Hilliard

A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation. 

## Compilation

Compilation requires CMake. Follow these steps to compile the code:

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

At the end of the process the executables `ch_1D`, `ch_2D`, and `ch_3D` will be placed in the `build/bin` folder.

## Usage

The executables run simulations in the corresponding spatial dimension and take one mandatory argument: a `TOML` file containing the simulation options.

Look at the `examples` folder for some runnable input files. Most of the options should be self-explanatory. The code prints an output consisting of four columns (time, free energy per bin, mass per bin and time step) both to the standard output and to file. It also appends configurations to the trajectory files (one for each species). In addition, the code also prints standalone configuration files every time 

The sections below describe the input keys. Non-mandatory options show their default values in brackets. TOML numeric values may generally be written as integers or floating-point values; arrays must be homogeneous.

- `steps` (integer, required): number of integration steps to run. Must be >= 0. Parsed as a 64-bit integer.
- `seed` (integer, optional, default: `time(NULL)`): seed used to initialise the RNG (parsed as 64-bit integer).
- `free_energy` (string, required): selects one of `landau`, `simple_wertheim`, `saleh`, `generic_wertheim`, or `ricci` (see [Free energy models](#free-energy-models)).
- `N` (integer, required): linear size per dimension. Must be a power of two. The total number of cells is `N` (1D), `N*N` (2D) or `N*N*N` (3D) depending on the executable.
- `k` (single value or array, required): interfacial penalty coefficient(s). The code accepts either a single numeric value (applied to all species) or an array with one value per species. Values may only be specified as floating point numbers.
- `dt` (double, required): time step used by the integrator.
- `dx` (double, optional, default: `1.0`): physical bin size (units used by the user). The code rescales `dx` internally according to `distance_scaling_factor`.
- `distance_scaling_factor` (double, optional, default: `1.0`): rescales user lengths to the internal units. Internally the code multiplies `dx` by `user_to_internal` and rescales `k` and densities accordingly; changing this can improve numerical stability for particular models.
- `integrator` (string, optional, default: `"euler"`): integration scheme. Supported values are `euler`, `euler_mobility`, `pseudospectral`, `pseudospectral_mobility`, and `bailo`. CUDA implementations are selected automatically where available.
- `use_CUDA` (bool, optional, default: `false`): enable CUDA-enabled integrators (when built with CUDA support).
- `species_evolution` (string or array, optional, default: `"CH"`): evolution equation for each species. Accepted values are `"CH"` (conserved Cahn–Hilliard dynamics) and `"AC"` (non-conserved Allen–Cahn dynamics). A single value is replicated for all species; an array supplies one value per species. Allen–Cahn is currently supported only by the CPU `euler` integrator.
- `species_AC_chemical_potential` (number or array, optional, default: `0.0`): per-species constant used by Allen–Cahn evolution. For an `AC` species, the update is `∂rho/∂t = -M (mu - species_AC_chemical_potential * rho)`. It is ignored for `CH` species.
- `mobility.type` (string, optional, default: `"constant"`): mobility implementation. Supported values are `constant`, `free_energy`, `regularised`, and `gel`; see [Mobility](#mobility).
- `mobility.M` (double, optional, default: `1.0`): base mobility used by `constant`, `free_energy`, and `regularised` mobility.
- `print_every` (integer, optional, default: `0`): frequency (in steps) at which the main energy/mass/time output line is written to `energy.dat` and to stdout. When `0` no periodic energy output is produced.
- `output.print_pressure` (bool, optional, default: `false`): if `true`, the average pressure is appended to the energy output and pressure fields are written to trajectory output.
- `output.print_pressure_strategy` (string, optional, default: `"linear"`): controls pressure trajectory printing. Supported values:
  - `"linear"`: print pressure files every `output.print_pressure_every` steps.
  - `"log"`: print pressure files at times round(`pressure_log_n0 * pressure_log_fact^N`) where `N` is the number of pressure frames already printed.
- `output.print_pressure_every` (integer, required for linear pressure output): frequency (in steps) to compute and write the pressure.
- `output.print_chemical_potential` (bool, optional, default: `false`): if `true`, the average chemical potential is appended to the energy output and chemical-potential fields are written to trajectory output.
- `output.print_chemical_potential_strategy` (string, optional, default: `"linear"`): controls chemical potential trajectory printing. Supported values:
  - `"linear"`: print chemical-potential files every `output.print_chemical_potential_every` steps.
  - `"log"`: print chemical potential files at times round(`chemical_potential_log_n0 * chemical_potential_log_fact^N`) where `N` is the number of chemical-potential frames already printed.
- `output.print_chemical_potential_every` (integer, required for linear chemical-potential output): frequency (in steps) to compute and write the chemical potential.
- `print_trajectory_strategy` (string, optional, default: `"linear"`): controls trajectory printing. Supported values:
  - `"linear"`: print configurations at fixed intervals using `print_trajectory_every`.
  - `"log"`: print configurations at times round(`log_n0 * log_fact^N`) where `N` is the number of trajectory frames already printed.
- `print_trajectory_every` (integer, optional, default: `0`): when using the `linear` strategy, append configurations to the trajectory every this many steps. If `0` no trajectory is appended.
- `print_last_every` (integer, optional): frequency (in steps) to write the `last_*` snapshot files. Defaults to the value of `print_trajectory_every` for the `linear` strategy. When using the `log` strategy `print_last_every` is required and must be explicitly provided.
- `output.log_n0` (integer, required for logarithmic trajectory output): base step for the logarithmic spacing.
- `output.log_fact` (double, required for logarithmic trajectory output): multiplicative factor for the logarithmic spacing.
- `output.pressure_log_n0` (integer, required for logarithmic pressure output): base step for pressure output spacing.
- `output.pressure_log_fact` (double, required for logarithmic pressure output): multiplicative factor for pressure output spacing.
- `output.chemical_potential_log_n0` (integer, required for logarithmic chemical-potential output): base step for chemical-potential output spacing.
- `output.chemical_potential_log_fact` (double, required for logarithmic chemical-potential output): multiplicative factor for chemical-potential output spacing.
- `output_path` (string, optional, default: `.`): directory where `last_*`, `init_*`, `*_*.dat`, and `energy.dat` are written.
- `output.path` (string, optional, default: `.`): directory used by the `Printer` for snapshots and field output. In normal runs keep this equal to `output_path`.
- `output.print_vtk` (bool, optional, default: `false`): when `true` produce VTK files instead of the native text format for snapshots.
- `output.trajectory_path` (string, optional): directory where trajectory files are written. When `output.print_vtk = true` this key is mandatory; otherwise it defaults to `output.path`.
- `load_from` (string, optional): path to a plain-text file used to initialise the fields. If present the file is parsed and used as the starting configuration (see "Initial configuration" below). When restarting from a `load_from` file the program will append to existing outputs.
- `initial_density` (single value or array, optional if `load_from` is present): average density used to randomly generate the initial configuration when `load_from` is not given. Accepts a single numeric value (applied to all species) or an array with one value per species.
- `initial_A` (double, optional, default: `1e-2`): amplitude of an optional sinusoidal modulation applied to the initial condition.
- `initial_N_peaks` (integer, optional, default: `0`): number of peaks for the initial sinusoidal modulation. When `0` the initial condition is purely random (white-noise like) around `initial_density`.

Notes about numeric fields and arrays:
- Keys parsed as per-species arrays may be provided either as a single value or as an array. A single value is replicated to the number of species when needed.
- `k`, `initial_density`, `species_evolution`, and `species_AC_chemical_potential` therefore support either one value for all species or one value per species.


## Initial configuration

The program accepts an initial-configuration file in the native text format produced by the program itself (or a file using the same layout). Alternatively, when `load_from` is not present the initial configuration is generated randomly from the TOML options described below.

Native text format (accepted by `load_from`)
- The file is plain text and may start with a header line in the form printed by the program:

  ```
  # step = <step>, t = <time>, size = Nx[ xNy[ xNz]]
  ```

- After the header the field values follow as whitespace-separated numbers. The layout depends on dimensionality (`N`):
  - 1D: the file contains `N` non-comment lines for each species; each line is a single number (the density for that bin).
  - 2D: the file contains `N` non-comment lines for each species; each line contains `N` whitespace-separated numbers (rows of the 2D grid, left-to-right).
  - 3D: the file contains `N*N` non-comment lines for each species; each line contains `N` whitespace-separated numbers. The data are written in a flattened order so that each line holds `N` values and the set of `N*N` lines represents the full `N x N x N` grid.

- Lines beginning with `#` are treated as comments and ignored when the file is read back with `load_from`.

Notes about `load_from` parsing
- The program reads the file species-by-species: for each species it skips comment lines and reads the expected number of data lines (see dims above). The native output produced by a run of the program (the `last_*` or `traj_*.dat` native files) can be reused as a `load_from` file for a later run. For multi-species systems, an initial configuration can be created by concatenating one file for each species, *e.g.* `cat last_?.dat > initial.dat`.

Generating the initial configuration randomly from the TOML input
- If `load_from` is not given, the initial field is created using the following TOML options:
  - `initial_density`: a single numeric value (applied to all species) or an array with one value per species. The code uses `_config_array_values<double>` to read this key and will replicate a single value to all species if needed.
  - `initial_A` (double, default `1e-2`): amplitude of an optional sinusoidal modulation applied across bins.
  - `initial_N_peaks` (integer, default `0`): number of peaks of the sinusoidal modulation. When `0` the initialization is purely random (white-noise like) around `initial_density`.

- The initialization algorithm (as implemented in the code) computes a modulation wavevector

  initial_k = 2 * pi * initial_N_peaks / N

  then for each linear bin index `bin` (from 0 to `grid_size` - 1) it computes

  modulation = initial_A * cos(initial_k * bin)

  and a `random_factor`:
  - if `initial_N_peaks == 0`: `random_factor = drand48() - 0.5` (zero-mean white noise)
  - otherwise: `random_factor = 1.0 + 0.02 * (drand48() - 0.5)` (small random perturbation around 1)

  For each species `i` the initial value placed in the internal grid is

  - if `average_rho != 0`: rho = average_rho * (1.0 + 2.0 * modulation * random_factor)
  - else: rho = 2.0 * modulation * random_factor

- The values are generated in user units and are later rescaled internally according to `distance_scaling_factor` (see the `distance_scaling_factor` / `dx` discussion above). The RNG seed can be set via `seed` (defaults to `time(NULL)`).

Practical tips
- If you want to produce a reproducible starting file, set the `seed` option or run the program once to write a `last_*` snapshot (or `traj_*.dat`), then reuse the output as `load_from` for subsequent runs.
- Make sure the file you provide as `load_from` matches the dimensionality (`ch_1D` / `ch_2D` / `ch_3D`) and the declared `N` used in the TOML configuration.

## Integrators

The code implements several time integrators. Choose one with the `integrator` key. All schemes use periodic boundary conditions. The following table summarizes the current support:

| Integrator | Non-constant mobility | Allen-Cahn | CUDA |
| --- | --- | --- | --- |
| `euler` | no | CPU only | yes |
| `euler_mobility` | yes | no | yes |
| `pseudospectral` | no | no | yes |
| `pseudospectral_mobility` | yes | no | no |
| `bailo` | no | no | no |

When `use_CUDA = true`, a CUDA implementation is selected where one exists. CUDA support must be enabled at build time.

- `euler` (default) — explicit finite-difference Euler stepping
  - Description: explicit finite-difference stepping. It supports both `CH` and `AC` species and computes local laplacians on the grid.
  - Use when: you want the simplest implementation or need Allen-Cahn dynamics. Explicit schemes require a sufficiently small `dt` for stability.

- `euler_mobility` — explicit Euler with variable mobility
  - Description: explicit finite-volume flux stepping with spatially varying mobility.
  - `mobility.with_noise` (bool, optional, default: `false`): add stochastic flux noise.
  - `mobility.noise_rescale_factor` (double, optional, default: `1.0`): multiply the noise amplitude by this factor.
  - Use when: mobility varies in space or stochastic fluxes are required. This integrator supports `CH` only.

- `pseudospectral` — semi-implicit pseudospectral (FFT) integrator with constant mobility
  - Description: FFT-based semi-implicit scheme for conserved dynamics.
  - Use when: you need better stability or larger time steps and can use FFTW on CPU or cuFFT on CUDA. This integrator supports `CH` only.

- `pseudospectral_mobility` — pseudospectral scheme with explicit treatment of mobility corrections
  - Description: semi-implicit FFT scheme with variable mobility. The mobility remainder is treated explicitly; the optional GMRES path solves the fully implicit variable-mobility correction.
  - `pseudospectral.use_gmres` (bool, optional, default: `false`): use the fully implicit variable-mobility correction.
  - `pseudospectral.gmres_restart` (integer, optional, default: `30`): GMRES restart length.
  - `pseudospectral.gmres_max_iter` (integer, optional, default: `200`): maximum GMRES iterations.
  - `pseudospectral.gmres_tol` (double, optional, default: `1e-10`): GMRES convergence tolerance.
  - Use when: mobility varies and the stability benefits of a semi-implicit spectral solver are useful. This integrator supports `CH` only and is CPU-only.

- `bailo` — implicit finite-volume integrator
  - Description: finite-volume scheme using expansive/contractive free-energy splitting and a nonlinear solve.
  - Use when: the selected free-energy model implements the required expansive and contractive derivatives. This integrator supports constant mobility and `CH` only.

### Allen-Cahn integration

Set `species_evolution = "AC"` for a one-species model, or provide one value per species for mixtures. With the CPU `euler` integrator, an `AC` species evolves according to

$$
\frac{\partial \rho}{\partial t} = -M\left(\mu - h\rho\right),
$$

where `h` is the corresponding value in `species_AC_chemical_potential`. Other species can remain conserved in the same run by using, for example, `species_evolution = ["CH", "AC"]`. Allen-Cahn is not currently available in the CUDA or pseudospectral implementations.

### S splitting parameter for pseudospectral integrators
- Key: `pseudospectral.S` (double, optional, default `0.0`)
- Key: `pseudospectral.use_dealias` (bool, optional, default `false`): apply the two-thirds spectral dealiasing filter.
- Role: `S` is a linear splitting parameter used in the semi-implicit spectral update. The pseudospectral integrator advances the Fourier components using a denominator of the form

  denom = 1 + dt * M * (S * k^2 + 2 * k_laplacian * k^4)

  and subtracts an implicit linear contribution `S * rho_hat` from the explicit free-energy derivative term. In practice this moves a linearised part of the chemical-potential derivative into the implicit side, increasing numerical stability for stiff nonlinearities.

- Guidance:
  - Default `S = 0.0` leaves the method semi-implicit only for the highest-order ($k^4$) term. Increasing `S` increases implicit stabilization and typically allows larger `dt` at similar stability.
  - Choose `S` based on the scale of the (local) derivative of the free-energy: a reasonable heuristic is to set `S` approximately equal to the largest expected value of $f''(\rho)$ (the derivative of $f$ with respect to density) in user units. For example, for a Landau model with $f' = -\epsilon * \psi + \psi^3$, $f'' = -\epsilon + 3 \psi^2$, since $|\psi| \approx 1$, for small values of $\epsilon$, $S \simeq 3$ is a sensible starting point.
  - Larger `S` stabilises the scheme but may overdamp fast modes and reduce accuracy; tune `S` and `dt` together.

Example (TOML):

```
[pseudospectral]
S = 2.0
use_dealias = true
```

## Mobility

Mobility is selected with `mobility.type` and is used by the integrator to construct the mobility field.

- `constant`: spatially uniform mobility `M`, supported by every integrator.
- `free_energy`: obtains mobility from the free-energy model's mobility implementation, with `M` as the model scale. Use a model that implements mobility coupling.
- `regularised`: regularises the mobility near zero density. Requires `mobility.rho_min`.
- `gel`: Landau gel mobility. Requires `mobility.phi_critical`, `mobility.c_0`, `mobility.M_c`, and `landau.epsilon`.

Non-constant mobility requires `euler_mobility` or `pseudospectral_mobility`.

## Free energy models

The code supports five free-energy models. Select one with the `free_energy` key: `landau`, `simple_wertheim`, `saleh`, `generic_wertheim`, or `ricci`.

### Landau free energy

Set `[landau].epsilon` or `[landau].T`; the model computes `epsilon = 1 - T` when `T` is supplied. `distance_scaling_factor` must remain `1.0` for this model.

This is the classic Ginzburg–Landau mean-field expression for the bulk free energy density:

$$
f_{\rm bulk}(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4
$$

where $\epsilon = (T_c - T) / T_c$ and $\psi$ is the order parameter of the phase transition. The total free-energy density used in the code is

$$
f(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4 + \kappa |\nabla \psi|^2
$$

where $\kappa$ defaults to one.

### Wertheim free energy

The model selected by `free_energy = "simple_wertheim"` requires `[wertheim].valence`, `[wertheim].B2`, and `[wertheim].delta`. `[wertheim].regularisation_delta` is optional and defaults to `0.0`.

This is the expression derived by Wertheim through his Thermodynamic Perturbation Theory to describe the thermodynamics of valence-limited fluids. The free energy density that is implemented in this code reads:

$$
\beta f_{\rm bulk} = \beta f_{\rm ref} + \beta f_{\rm bond}
$$

where $\beta f_{\rm ref} = \rho \ln(\rho) - \rho + B_2 \rho^2$, with $B_2$ second virial coefficient, is the free energy of the reference system (*i.e.* the system where no bonding is possible) and 

$$
\beta f_{\rm bond} = M \rho \left(\ln(X) + \frac{1}{2} (1 - X) \right)
$$

is the free energy that takes into account bonding. Here $M$ is the valence of each particle and $X$ is the probability that a patch is unbound. The latter can be estimated through a law of mass-action and is equal to

$$
X(\rho) = \frac{-1 + \sqrt{1 + 4 M \Delta \rho}}{2 M \Delta \rho}
$$

where $\Delta = v_b e^{-\Delta G / R T}$, $v_b = 1.6606$ nm $^3$ and $\Delta G$ is the DNA hybridisation free energy.

### Saleh free energy

Select this model with `free_energy = "saleh"`. It requires `[saleh].B2`, `[saleh].delta_AA`, and `[saleh].delta_BB`. `[saleh].B3` defaults to `0.0`, and `[saleh].valence` is a scalar or three-element array defaulting to `3`.

This is the Wertheim free energy for a ternary mixture of valence-limited particles. The three species $A$, $B$ and $C$ have the same intra- and inter-species repulsion (provided by a second virial coefficient that takes the same value for every interaction). However, $A$ can bind only to $A$ or to half of the sites on $C$, $B$ only to $B$ or to half of the sites of $C$, so that the $C$ species acts as a linker.

### Ricci Wertheim free energy

Select this model with `free_energy = "ricci"`. It requires `[ricci].B2`, `[ricci].delta_00`, and `[ricci].delta_12`.

### Generic Wertheim free energy

Select this model with `free_energy = "generic_wertheim"`. It requires:

- one `[[generic_wertheim.species]]` table per species, each with a `patches` integer array;
- `[[generic_wertheim.deltas]]` tables with an `interaction` such as `"0-1"` and parameters accepted by the `Delta` parser (`T`, `deltaH`, `deltaS`, `salt`, `sticky_size`, or a direct `value`);
- `[[generic_wertheim.B2s]]` tables with an `interaction` and numeric `value`.

The number of `B2s` entries is normally $N(N+1)/2$. Set `generic_wertheim.allow_unspecified_B2s = true` to run with an incomplete set; unspecified coefficients remain zero.
