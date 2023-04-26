# Cahn-Hilliard

A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation. 

## Landau free energy

By default the code uses the classic Ginzburg–Landau mean-field expression for the bulk free energy density:

$$
f_{\rm bulk}(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4
$$

where $\epsilon = (T_c - T) / T_c$ and $\psi$ is the order parameter of the phase transition. The total free-energy density used in the code is

$$
f(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4 + \kappa |\nabla \psi|^2
$$

where $\kappa$ defaults to one.

## Wertheim free energy

By using the `--free-energy wertheim` switch (see [usage](#usage)) it is possible to use the free energy of a system made of tetrafunctional DNA arms (with parameters taken from [here](https://doi.org/10.1021/acsnano.6b08287)) modelled as tetrafunctional patchy particles through Wertheim's theory. Here we use the following free energy density:

$$
\beta f_{\rm bulk} = \beta f_{\rm ref} + \beta f_{\rm bond}
$$

where $\beta f_{\rm ref} = \rho \ln(\rho) - \rho + B_2 \rho^2$  with $B_2 = 2100$ nm $^{3}$ second virial coefficient, is the free energy of the reference system (*i.e.* the system where no bonding is possible) and 

$$
\beta f_{\rm bond} = M \rho \left(\ln(X) + \frac{1}{2} (1 - X) \right)
$$

is the free energy that takes into account bonding. Here $M = 4$ is the valence of each particle and $X$ is the probability that a patch is unbound. The latter can be estimated through a law of mass-action and is equal to

$$
X(\rho) = \frac{-1 + \sqrt{1 + 4 M \Delta \rho}}{2 M \Delta \rho}
$$

where $\Delta = v_b e^{-\Delta G / R T}$, $v_b = 1.6606$ nm $^3$ and $\Delta G$ is the DNA hybridisation free energy, here estimated through Santa Lucia's nearest-neighbour model.

## Compilation

Compilation requires CMake. Follow these steps to compile the code:

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

At the end of the process an output executable, `ch`, will be placed in the `build` folder.

## Usage

```
  cahn-hilliard [OPTION...]

  -s, --steps arg        Number of iterations
  -N arg                 The size of the square grid (default: 64)
  -f, --free-energy arg  The bulk free energy expression to be used (supported values are 'landau' and 'wertheim') (default: landau)
  -e, --epsilon arg      The distance from the critical point in the 'landau' free energy (default: 0.9)
  -T, --temperature arg  Temperature (in Kelvin), used by the 'wertheim' free energy (default: 300)
  --dt arg               The integration time step (default: 0.01)
  -M arg                 The transport coefficient M of the Cahn-Hilliard equation (default: 1.0)
  -H arg                 The size of the mesh cells (default: 1.0)
  -a, --average-psi arg  Average value of the order parameter (default: 0)
  --psi-noise arg        Random noise for the initial psi (default: 0.1)
  -p, --print-every arg  Number of iterations every which the state of the system will be appended to the trajectory.dat file (0 means never) (default: 0)
  -k arg                 Strength of the interfacial term of the Cahn-Hilliard equation (default: 1.0)
  -h, --help             Print usage
```
