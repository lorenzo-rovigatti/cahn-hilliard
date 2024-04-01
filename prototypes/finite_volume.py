# some of the code has been taken from https://github.com/sergiopperez/Cahn_Hilliard_Finite_Volume

import numpy as np
from scipy import optimize

B2 = 10600
valence = 4
Delta = 59769.7978988913
two_valence_delta = 2 * valence * Delta

def F_der_contractive(rho):
    F_der = np.log(rho) + 2.0 * B2 * rho

    return F_der

def F_der_expansive(rho):
    X = (-1.0 + np.sqrt(1.0 + 2.0 * two_valence_delta * rho)) / (two_valence_delta * rho)

    F_der = -valence * np.log(X)

    return F_der

def to_solve(rho_curr, rho_next, epsilon, dx, dt, M, n):
    F_der_con = F_der_contractive(rho_next)
    F_der_exp = F_der_expansive(rho_curr)

    Lap = epsilon**2 * (np.roll(rho_curr, 1) - 2 * rho_curr + np.roll(rho_curr, -1) + np.roll(rho_next, 1) - 2 * rho_next + np.roll(rho_next, -1)) / (2 * dx**2)
    Fhalf = -(np.roll(F_der_con - F_der_exp - Lap, 1) - F_der_con + F_der_exp + Lap) / dx

    return rho_next - rho_curr + (Fhalf - np.roll(Fhalf, -1)) * dt / dx


def main():
    rho0 = np.fromfile('init_0.dat', sep=" ")
    n = rho0.shape[0] # Number of cells
    dx = 20 # Width of cells
    x = np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
    epsilon = np.sqrt(2 * 1e7) # Parameter epsilon   
    dt = 1e-3 # Time step
    M = 1
    ntimes = 10000
    
    rho = np.zeros([n,ntimes+1]) # Density matrix
    rho[:,0] = rho0 # First column of density matrix is initial density

    for i in np.arange(ntimes): #Temporal loop
        
        rho[:, i+1], infodict, ier, mesg = optimize.fsolve(lambda rhon: to_solve(rho[:,i], rhon, epsilon, dx, dt, M, n), rho[:,i], full_output = True)
        rho_curr = rho[:, i + 1]

        if i % 10 == 0:
            X = (-1.0 + np.sqrt(1.0 + 2.0 * two_valence_delta * rho_curr)) / (two_valence_delta * rho_curr)
            grad = (np.roll(rho_curr, -1) - rho_curr) / dx
            F_bulk = rho_curr * (np.log(rho_curr) - 1.0) + B2 * rho_curr**2 + valence * rho_curr * (np.log(X) + 0.5 * (1.0 - X))
            F_interf = epsilon**2 / 2 * grad**2 # works in 1D
            F = np.sum(F_bulk + F_interf) * dx**3 / n

            print(f"{dt * i:.2f} {F:.5f}")
    
    with open("last_rho.dat", "w") as f:
        rho[:,-1].tofile(f, "\n")
        f.write("\n")
    

if __name__ == '__main__':
    main()