# Cahn-Hilliard Simulation Architecture

## Overview

This codebase implements numerical simulations of spinodal decomposition through the Cahn-Hilliard equation. The architecture is built around a **template-based dimensional abstraction** (`1D`, `2D`) and a **plugin-based design** for swappable physics models, integration schemes, and mobility functions.

## Project Root Structure

```
src/                  # Main simulation code
├── main.cpp          # Entry point with Manager class
├── CahnHilliard.h/cpp # Core PDE solver and grid management
├── SimulationState.h  # Shared simulation state across components
├── Object.h/cpp       # Base class for all managed objects (logging, config parsing)
├── models/            # Free energy models (pluggable)
├── integrators/       # Time stepping schemes (pluggable)
├── mobility/          # Mobility models (pluggable)
└── utils/             # Grid representations and utilities
```

## Core Class Hierarchy

### 1. Base Class: `Object`

- **Location**: [src/Object.h](src/Object.h)
- **Purpose**: Base class for all configurable, loggable components
- **Key Features**:
  - Provides TOML configuration parsing utilities (`_config_value<T>`, `_config_optional_value<T>`, etc.)
  - Integrates with spdlog for structured logging (info, warning, critical)
  - Enables polymorphic management of components

All major classes inherit from `Object`:
- `Manager`
- `CahnHilliard<dims>`
- `Integrator<dims>`
- `FreeEnergyModel` (and its subclasses)
- `IMobility` (and its subclasses)

### 2. Manager: Orchestration Layer

- **Location**: [src/main.cpp](src/main.cpp) - nested class in `main()`
- **Responsibilities**:
  - **Configuration Loading**: Parses TOML file, validates output directory
  - **Component Initialization**: 
    - Instantiates the appropriate free energy model based on config
    - Creates the `CahnHilliard<DIM>` solver
    - Sets up trajectory output files
  - **Simulation Loop**: 
    - Manages time-stepping across `_steps` iterations
    - Coordinates printing of mass, energy, pressure, and trajectory data
    - Handles both linear and logarithmic printing strategies
  - **Restart Support**: Can resume from saved configurations via `load_from` parameter

**Key Members**:
- `_system`: Unique pointer to the `CahnHilliard<DIM>` solver
- `_sim_state`: Central simulation state (shared with all components)
- `_trajectories`: Output streams for species density trajectories

### 3. SimulationState: Data Container

- **Location**: [src/SimulationState.h](src/SimulationState.h)
- **Purpose**: Centralized state shared across all components
- **Key Members**:
  - `rho`: Density field for all species (MultiField<double>)
  - `mobility`: Mobility field (may vary in space/time depending on model)
  - `model`: Unique pointer to the current free energy model
  - `CUDA_rho`, `CUDA_mobility`: GPU device pointers (when using CUDA)
  - `time_step`: Current simulation time step
  - `use_CUDA`: Flag enabling GPU computation

All components receive a reference to `SimulationState` and coordinate through it.

### 4. CahnHilliard<dims>: PDE Solver

- **Location**: [src/CahnHilliard.h](src/CahnHilliard.h)
- **Template Parameter**: `dims` (1 or 2 for 1D and 2D simulations)
- **Purpose**: High-level solver managing the grid and delegating time integration

**Key Responsibilities**:
- Grid discretization: Converts between coordinates and linear indices
- Gradient computation in real space (finite differences)
- Laplacian computation
- Output formatting: Writes species densities, pressure, energy to files/streams
- Physics calculation: Computes average mass, free energy, pressure

**Key Methods**:
- `evolve()`: Advances one time step by calling the integrator
- `average_mass()`, `average_free_energy()`, `average_pressure()`: Compute statistics
- `print_*()`: Output to files or streams

**Key Members**:
- `model`: Pointer to free energy model
- `integrator`: Pointer to the time stepping scheme (swappable)
- `mobility`: Unique pointer to mobility model (swappable)
- Various grid parameters: `N`, `grid_size`, `dx`, `dt`, `k_laplacian`

### 5. Free Energy Models (Pluggable)

- **Location**: [src/models/](src/models/)
- **Base Class**: `FreeEnergyModel`

**Abstract Interface**:
- `N_species()`: Number of chemical species
- `der_bulk_free_energy()`: Compute derivative of free energy density (CPU version)
- `der_bulk_free_energy_*()`: Expansive/contractive splitting variants
- `pressure()`: Optional pressure calculation (for pressure-coupled models)
- `set_mobility()`: Optional model-specific mobility (e.g., concentration-dependent)

**Concrete Implementations**:
1. **`Landau`** ([src/models/Landau.h](src/models/Landau.h)): Simple one-component model
   - $f(\rho) = \frac{\epsilon}{2}(\rho - \rho_c)^2$
   
2. **`SimpleWertheim`** ([src/models/SimpleWertheim.h](src/models/SimpleWertheim.h)): Two-component associative liquid model
   
3. **`SalehWertheim`** ([src/models/SalehWertheim.h](src/models/SalehWertheim.h)): Enhanced Wertheim model
   
4. **`GenericWertheim`** ([src/models/GenericWertheim.h](src/models/GenericWertheim.h)): Flexible N-component Wertheim
   
5. **`RicciWertheim`** ([src/models/RicciWertheim.h](src/models/RicciWertheim.h)): Alternative Wertheim formulation

Each model encodes different thermodynamic behavior but provides the same interface for grid-based evaluation.

### 6. Integrators (Pluggable)

- **Location**: [src/integrators/](src/integrators/)
- **Base Class**: `Integrator<dims>` (template on spatial dimensions)

**Abstract Interface**:
- `evolve()`: Advance one time step
- `validate()`: Check configuration consistency
- `sync()`: Optional GPU ↔ CPU synchronization

**Concrete Implementations**:

1. **`EulerCPU<dims>`** ([src/integrators/EulerCPU.h](src/integrators/EulerCPU.h)): Explicit finite-difference Euler
   - Naive grid-based discretization of Cahn-Hilliard PDE
   - Limited by CFL condition on time step
   
2. **`EulerMobilityCPU<dims>`** ([src/integrators/EulerMobilityCPU.h](src/integrators/EulerMobilityCPU.h)): Euler with variable mobility
   - Extends EulerCPU to handle spatially/temporally varying mobility
   
3. **`PseudospectralCPU<dims>`** ([src/integrators/PseudospectralCPU.h](src/integrators/PseudospectralCPU.h)): FFT-based semi-implicit
   - Uses FFTW for efficient Laplacian computation in Fourier space
   - Semi-implicit time stepping for stability
   - Supports mobility splitting: $M = M_0 + (M - M_0)$
   - More efficient and stable than Euler, especially for long runs
   
4. **`PseudospectralMobilityCPU<dims>`** ([src/integrators/PseudospectralMobilityCPU.h](src/integrators/PseudospectralMobilityCPU.h)): Pseudospectral with variable mobility
   - Combines FFT Laplacian with finite-difference mobility splitting
   
5. **`BailoFiniteVolume<dims>`** ([src/integrators/BailoFiniteVolume.h](src/integrators/BailoFiniteVolume.h)): Finite volume scheme
   - Alternative discretization approach

Integrators are selected at runtime via the `integrator` config option.

### 7. Mobility Models (Pluggable)

- **Location**: [src/mobility/](src/mobility/)
- **Base Class**: `IMobility`

**Abstract Interface**:
- `update_mobility()`: Update the mobility field at current state

**Concrete Implementations**:

1. **`ConstantMobility`** ([src/mobility/ConstantMobility.h](src/mobility/ConstantMobility.h)): Spatially uniform
   - $M(\rho) = M_0$ (constant)
   
2. **`FreeEnergyMobility`** ([src/mobility/FreeEnergyMobility.h](src/mobility/FreeEnergyMobility.h)): Derived from free energy
   - Delegates to the model's `set_mobility()` method
   
3. **`RegularisedMobility`** ([src/mobility/RegularisedMobility.h](src/mobility/RegularisedMobility.h)): Regularized form
   
4. **`GelMobility`** ([src/mobility/GelMobility.h](src/mobility/GelMobility.h)): Specialized for gel systems

Mobility is selected at runtime via the `mobility.type` config option.

### 8. Utilities

- **Location**: [src/utils/](src/utils/)

Key components:

1. **`MultiField<T>`** ([src/utils/MultiField.h](src/utils/MultiField.h)): Generic multi-species grid representation
   - Stores multiple scalar fields (one per species) in a single buffer
   - Provides species views for convenient slicing
   - Supports both CPU and GPU memory layouts
   
2. **`SpeciesView<T>`** ([src/utils/MultiField.h](src/utils/MultiField.h)): Lightweight read-only view of a single species
   - Enables range-based for loops over species data
   
3. **`Gradient<dims>`** ([src/utils/Gradient.h](src/utils/Gradient.h)): Vector field storage
   - Stores directional derivatives for each point
   
4. **`Delta`** ([src/utils/Delta.h](src/utils/Delta.h)): Used to parse and store values for the delta parameter Wertheim's theory
   
5. **String utilities** ([src/utils/strings.h](src/utils/strings.h)): TOML string parsing and manipulation

## Data Flow: A Simulation Time Step

1. **Manager Loop** calls `_system->evolve()`
   
2. **CahnHilliard::evolve()**
   - Delegates to the selected `Integrator<dims>::evolve()`
   
3. **Integrator::evolve()** (e.g., PseudospectralCPU)
   - Reads current density from `SimulationState::rho`
   - Calls `FreeEnergyModel::der_bulk_free_energy()` to get chemical potential
   - Applies mobility from `IMobility::update_mobility()` and `SimulationState::mobility`
   - Computes Laplacian (via FFT in Fourier space)
   - Advances `SimulationState::rho` to next time step
   
4. **Optionally**: Output statistics and fields to files

## Configuration (TOML)

All components receive configuration from a TOML file. Key sections:

```toml
[general]
steps = 10000
N = 512
dt = 0.001
dx = 1.0
free_energy = "landau"           # or "simple_wertheim", "generic_wertheim", etc.
integrator = "pseudospectral"    # or "euler", "finite_volume", etc.

[mobility]
type = "constant"                # or "free_energy", "regularised", etc.
M0 = 1.0

[output]
print_every = 100
print_trajectory_every = 1000
output_path = "."

[model_specific]
# Varies by free_energy choice (e.g., epsilon for Landau)
epsilon = 0.001
```

The `Object` base class provides utilities for type-safe configuration access with optional defaults.

## Compilation & Build

- **Build System**: CMake
- **Conditional Features**:
  - `NOCUDA` flag disables GPU code paths
  - Separate build directories for CPU (`build/`) and CUDA (`cuda/`)
- **Dependencies**:
  - TOML++ for configuration parsing
  - FFTW3 for FFT-based integrators
  - spdlog for logging
  - Cxxopts for command-line parsing (examples)

## Key Design Patterns

1. **Strategy Pattern**: Integrators, mobility models, and free energy models are swappable at runtime
2. **Template Method**: Integrator base class defines the validation flow; subclasses implement `evolve()`
3. **Visitor Pattern**: MultiField provides SpeciesView for iteration over species
4. **Dependency Injection**: Manager creates components and passes references to SimulationState
5. **RAII**: Unique pointers manage object lifetimes; file streams are automatically closed

## Extension Points

To add new functionality:

1. **New Free Energy Model**: Inherit from `FreeEnergyModel`, implement `der_bulk_free_energy()` and `N_species()`
   - Register in `Manager` constructor
   
2. **New Integrator**: Inherit from `Integrator<dims>`, implement `evolve()`
   - Select via `integrator` config option
   
3. **New Mobility**: Inherit from `IMobility`, implement `update_mobility()`
   - Select via `mobility.type` config option
   
4. **New Output Statistics**: Add method to `CahnHilliard` and call from `Manager::run()`

## Performance Considerations

- **FFT Integrators**: Preferred for long simulations; requires periodic boundary conditions
- **Finite Difference**: More flexible geometry; limited by explicit time stepping stability
- **GPU Acceleration**: Available for free energy derivatives (CUDA code in `src/CUDA/`)
- **Grid Dimensions**: 1D is fastest; 2D requires $N^2$ memory; support up to 2D currently
