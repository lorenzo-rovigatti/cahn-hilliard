# Cahn-Hilliard

A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation. The code uses the classic Ginzburg–Landau mean-field expression for the bulk free energy density:

$$
f_{\rm bulk}(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4
$$

where $\epsilon = (T_c - T) / T_c$ and $\psi$ is the order parameter of the phase transition. The total free-energy density used in the code is thus

$$
f(\psi) = -\frac{1}{2}\epsilon \psi^2 + \frac{1}{4} \psi^4 + \frac{1}{2} |\nabla \psi|^2
$$

## Compilation and usage

Use `compile.sh`, which assumes that g++ is installed and in the PATH, to compile the code, and run it with

```
$ ./ch <N> <epsilon> <psi_average> <steps>
```

where

* `N` is the size of the square grid that will be used to discretise the problem. This value should be a power of two.
* `epsilon` is the parameter controlling the distance from the critical point (see above). Phase separation requires $\epsilon > 0$. The larger the number, the farther from the critical point.
* `psi_average` controls the initial conditions of the grid. Each cell of the grid will be assigned a random value that is `psi_average` $\pm 0.5$.
* `steps` is the number of iterations that the simulation will run for.

At the end of the simulation the final grid will be printed to the file `last.dat`.
