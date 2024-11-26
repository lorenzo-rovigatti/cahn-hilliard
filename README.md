# Cahn-Hilliard

A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation. 

## Compilation

Compilation requires CMake. Follow these steps to compile the code:

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

At the end of the process two executables, `ch_1D` and `ch_2D`, will be placed in the `build` folder.

## Usage

The two executables are used to run simulations in 1D and 2D, respectively, and take a single mandatory option, which is a `TOML` file containing the options specifying the behaviour of the simulation.

Look at the `examples` folder for some runnable input files. Most of the options should be self-explanatory.

## Free energy models

The code supports three free energy models, which can be specified by setting the `free_energy` key in the input to `landau`, `simple_wertheim` or `saleh`.

## Landau free energy

This is the classic Ginzburg–Landau mean-field expression for the bulk free energy density:

$$
f_{\rm bulk}(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4
$$

where $\epsilon = (T_c - T) / T_c$ and $\psi$ is the order parameter of the phase transition. The total free-energy density used in the code is

$$
f(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4 + \kappa |\nabla \psi|^2
$$

where $\kappa$ defaults to one.

## Wertheim free energy

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

## Saleh free energy

This is the Wertheim free energy for a ternary mixture of valence-limited particles. The three species $A$, $B$ and $C$ have the same intra- and inter-species repulsion (provided by a second virial coefficient that takes the same value for every interaction). However, $A$ can bind only to $A$ or to half of the sites on $C$, $B$ only to $B$ or to half of the sites of $C$, so that the $C$ species acts as a linker.
