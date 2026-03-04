# Cahn-Hilliard

A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation. 

## Compilation

Compilation requires CMake. Follow these steps to compile the code:

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

At the end of the process three executables, `ch_1D`, `ch_2D`, and `ch_3D` will be placed in the `build` folder.

## Usage

The two executables are used to run simulations in 1D and 2D, respectively, and take a single mandatory option, which is a `TOML` file containing the options specifying the behaviour of the simulation.

Look at the `examples` folder for some runnable input files. Most of the options should be self-explanatory. The code prints an output consisting of four columns (time, free energy per bin, mass per bin and time step) both to the standard output and to file. It also appends configurations to the trajectory files (one for each species). In addition, the code also prints standalone configuration files every time 

Here is a (code-accurate) list of input keys and their behaviour. Non-mandatory options show their default values in brackets. Where a numeric value is accepted both as integer or floating point the code will accept either (integers are converted to the appropriate floating type when needed).

- `steps` (integer, required): number of integration steps to run. Must be >= 0. Parsed as a 64-bit integer.
- `print_every` (integer, optional, default: `0`): frequency (in steps) at which the main energy/mass/time output line is written to `energy.dat` and to stdout. When `0` no periodic energy output is produced.
- `print_average_pressure` (bool, optional, default: `false`): if `true`, the average pressure is appended to the energy output and written to `pressure.dat` when `print_pressure_every` > 0.
- `print_pressure_every` (integer, optional, default: `0`): frequency (in steps) to compute and write the pressure to `pressure.dat`.
- `print_trajectory_strategy` (string, optional, default: `"linear"`): controls trajectory printing. Supported values:
  - `"linear"`: print configurations at fixed intervals using `print_trajectory_every`.
  - `"log"`: print configurations at times round(`log_n0 * log_fact^N`) where `N` is the number of trajectory frames already printed.
- `print_trajectory_every` (integer, optional, default: `0`): when using the `linear` strategy, append configurations to the trajectory every this many steps. If `0` no trajectory is appended.
- `print_last_every` (integer, optional): frequency (in steps) to write the `last_*` snapshot files. Defaults to the value of `print_trajectory_every` for the `linear` strategy. When using the `log` strategy `print_last_every` is required and must be explicitly provided.
- `log_n0` (integer, required for `log` strategy): base step for the logarithmic spacing.
- `log_fact` (double, required for `log` strategy): multiplicative factor for the logarithmic spacing.
- `seed` (integer, optional, default: `time(NULL)`): seed used to initialise the RNG (parsed as 64-bit integer).
- `free_energy` (string, required): selects the free-energy model. Accepted values include `landau`, `simple_wertheim`, `saleh`, `generic_wertheim`, `ricci` (see the "Free energy models" section).
- `N` (integer, required): linear size per dimension. Must be a power of two. The total number of cells is `N` (1D), `N*N` (2D) or `N*N*N` (3D) depending on the executable.
- `k` (single value or array, required): interfacial penalty coefficient(s). The code accepts either a single numeric value (applied to all species) or an array with one value per species. Values may only be specified as floating point numbers.
- `dt` (double, required): time step used by the integrator.
- `dx` (double, optional, default: `1.0`): physical bin size (units used by the user). The code rescales `dx` internally according to `distance_scaling_factor`.
- `distance_scaling_factor` (double, optional, default: `1.0`): rescales user lengths to the internal units. Internally the code multiplies `dx` by `user_to_internal` and rescales `k` and densities accordingly; changing this can improve numerical stability for particular models.
- `integrator` (string, optional, default: `"euler"`): integration scheme. Supported values include `euler`, `euler_mobility`, `pseudospectral`, `pseudospectral_mobility`, `bailo`. When `use_CUDA = true` some integrators have CUDA implementations (the code chooses the appropriate variant automatically).
- `use_CUDA` (bool, optional, default: `false`): enable CUDA-enabled integrators (when built with CUDA support).
- `output.path` (string, optional, default: `.`): directory where `last_*`, `*_*.dat` and `energy.dat` are written.
- `output.print_vtk` (bool, optional, default: `false`): when `true` produce VTK files instead of the native text format for snapshots.
- `output.trajectory_path` (string, optional): directory where trajectory files are written. When `output.print_vtk = true` this key is mandatory; otherwise it defaults to `output.path`.
- `load_from` (string, optional): path to a plain-text file used to initialise the fields. If present the file is parsed and used as the starting configuration (see "Initial configuration" below). When restarting from a `load_from` file the program will append to existing outputs.
- `initial_density` (single value or array, optional if `load_from` is present): average density used to randomly generate the initial configuration when `load_from` is not given. Accepts a single numeric value (applied to all species) or an array with one value per species.
- `initial_A` (double, optional, default: `1e-2`): amplitude of an optional sinusoidal modulation applied to the initial condition.
- `initial_N_peaks` (integer, optional, default: `0`): number of peaks for the initial sinusoidal modulation. When `0` the initial condition is purely random (white-noise like) around `initial_density`.

Notes about numeric fields and arrays:
- Keys parsed with the helper `_config_array_values<T>` may be provided either as a single value or as an array. If a single value is provided and an `output_size` (for example the number of species) is required the single value is replicated to match the expected size.
- Arrays must be homogeneous (all elements of the same TOML numeric type).


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

The code implements several time integrators. Choose one with the `integrator` key in the TOML input. When `use_CUDA = true` the program will automatically select CUDA implementations if available.

- `euler` (default) — explicit finite-difference Euler stepping
  - Description: a straightforward explicit finite-difference discretisation of the continuity equation used for Cahn–Hilliard evolution. It computes local gradients and laplacians on the grid and advances densities explicitly by `dt`.
  - Use when: you want a simple, robust CPU implementation. Time-step `dt` must be chosen small enough for stability (explicit schemes are conditionally stable).

- `euler_mobility` — explicit Euler with (possibly variable) mobility
  - Description: same spatial discretisation as `euler` but supports non-constant, density-dependent mobility fields. The integrator optionally adds stochastic noise to the mobility-driven flux when configured.
  - Config/notes: this integrator advertises support for non-constant mobility in the code; it accepts mobility-related parameters through the model/mobility interfaces. Noise can be enabled via integrator-specific TOML keys (see model/integrator docs or inspect `EulerMobilityCPU` for exact keys).
  - Use when: your free-energy model requires spatially varying mobility or you want to include mobility noise.

- `pseudospectral` — semi-implicit pseudospectral (FFT) integrator with mobility splitting
  - Description: FFT-based semi-implicit scheme. The integrator advances the solution in Fourier space using a semi-implicit factor for the highest-order (Laplacian) terms, while lower-order and mobility-correction terms are handled explicitly in real space and added as a correction. This improves stability and allows larger `dt` than explicit Euler for the same spatial resolution.
  - Important config keys (examples):
    - `mobility.M0` (double, optional): base mobility used in the splitting (defaults to a value derived from the mobility object).
    - `semi_implicit.rho_floor` (double, optional, default `0.0`): clamp applied to `rho` before calling free-energy derivatives.
    - `semi_implicit.dealias` (bool, optional, default `false`): enable dealiasing in spectral operations.
  - Use when: you need better stability / larger time-steps and can pay the FFT cost. Works well with smooth fields and periodic boundary conditions.

- `pseudospectral_mobility` — pseudospectral scheme with explicit treatment of mobility corrections
  - Description: variant of the `pseudospectral` integrator that explicitly accounts for variable mobility via a splitting strategy (compute correction terms in real space and transform them when needed).
  - Use when: mobility varies significantly and you still want the stability benefits of a semi-implicit spectral solver.

### Notes and selection guidance
- Explicit schemes (`euler`, `euler_mobility`) are simple and sometimes faster per-step on small grids, but require smaller `dt` for stability.
- Pseudospectral integrators require FFT libraries (the code uses FFTW on CPU and cuFFT on CUDA). They permit larger `dt` and can be more efficient on large grids, especially when paired with optimized FFT backends.
- CUDA implementations exist for most integrators; enable them with `use_CUDA = true` and ensure the code was built with CUDA support.

### S splitting parameter for pseudospectral integrators
- Key: `pseudospectral.S` (double, optional, default `0.0`)
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

## Free energy models

The code supports three free energy models, which can be specified by setting the `free_energy` key in the input to `landau`, `simple_wertheim`, `saleh`, or `generic_wertheim`.

### Landau free energy

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

This is the Wertheim free energy for a ternary mixture of valence-limited particles. The three species $A$, $B$ and $C$ have the same intra- and inter-species repulsion (provided by a second virial coefficient that takes the same value for every interaction). However, $A$ can bind only to $A$ or to half of the sites on $C$, $B$ only to $B$ or to half of the sites of $C$, so that the $C$ species acts as a linker.
