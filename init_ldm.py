import sys
import numpy as np

if len(sys.argv) < 4:
    print("Usage is %s N tetramer-density R" % sys.argv[0], file=sys.stderr)
    exit(1)
    
N = int(sys.argv[1])
tetramer_rho = float(sys.argv[2])
R = float(sys.argv[3])
linker_rho = 2 * tetramer_rho / (1.0 + R);

A_density = 2 * tetramer_rho * (1.0 + np.random.normal(loc=1.0, size=(N, N)) * 1e-2)
A_density[0:int(N/2)] *= 0.01

B_density = 2 * tetramer_rho * (1.0 + np.random.normal(loc=1.0, size=(N, N)) * 1e-2)
B_density[int(N/2):N] *= 0.01

AA_density = linker_rho * (1.0 + np.random.normal(loc=1.0, size=(N, N)) * 1e-2);
BB_density = linker_rho * (1.0 + np.random.normal(loc=1.0, size=(N, N)) * 1e-2);

AB_density = 2 * R * linker_rho * (1.0 + np.random.normal(loc=1.0, size=(N, N)) * 1e-2);

with open("generated.dat", "w") as out:
    np.savetxt(out, A_density)
    np.savetxt(out, B_density)
    np.savetxt(out, AA_density)
    np.savetxt(out, BB_density)
    np.savetxt(out, AB_density)
    