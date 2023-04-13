# Cahn-Hilliard

A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation. The code uses the classic Ginzburg–Landau mean-field expression for the bulk free energy density:

$$
f_{\rm bulk}(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4
$$

where $\epsilon = (T_c - T) / T_c$ and $\psi$ is the order parameter of the phase transition. The total free-energy density used in the code is

$$
f(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4 + \frac{1}{2} \kappa |\nabla \psi|^2
$$

where $\kappa$ defaults to one.

## Compilation

Compilation requires CMake. Follow these steps to compile the code:

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

At the end of the process an output executable, `ch`, will be placed in the `build` folder.

## Usage

```
chan-hilliard [OPTION...]

-N arg                 The size of the square grid (default: 64)
-e, --epsilon arg      The distance from the critical point (default: 0.9)
    --dt arg           The integration time step (default: 0.01)
-a, --average-psi arg  Average value of the order parameter (default: 0)
-s, --steps arg        Number of iterations
-p, --print-every arg  Number of iterations every which the state of the system will be appended to the trajectory.dat file  (0 means never) (default: 0)
-k arg                 Strength of the interfacial term of the  Cahn-Hilliard equation (default: 1.0)
-h, --help             Print usage
```
