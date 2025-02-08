# some of the code has been taken from https://github.com/sergiopperez/Cahn_Hilliard_Finite_Volume

import numpy as np
from scipy.fftpack import fft, ifft
from scipy import optimize

B2 = 2190
valence = 4
Delta = 29696.93101147855
two_valence_delta = 2 * valence * Delta
rho_min = 1e-9
C = -1

def free_energy_density(rho):
    X = (-1.0 + np.sqrt(1.0 + 2.0 * two_valence_delta * rho)) / (two_valence_delta * rho)

    F_ref = rho * np.log(np.maximum(rho, rho_min)) - rho + B2 * rho**2
    F_bond = valence * rho * (np.log(X) + 0.5 * (1 - X))

    return F_ref + F_bond

def F_der_contractive(rho):
    F_der = np.log(rho) + 2.0 * B2 * rho

    return F_der

def F_der_expansive(rho):
    X = (-1.0 + np.sqrt(1.0 + 2.0 * two_valence_delta * rho)) / (two_valence_delta * rho)

    F_der = -valence * np.log(X)

    return F_der

def main():
    rho0 = np.loadtxt('init_0.dat', comments="#")
    n = rho0.shape[0] # Number of cells
    dx = 10 # Width of cells
    epsilon = np.sqrt(2 * 1e6) # Parameter epsilon   
    dt = 1e-6 # Time step
    M = 1
    ntimes = 1000
    
    rho = np.zeros([n,ntimes+1]) # Density matrix
    rho[:, 0] = rho0 # First column of density matrix is initial density

    # Wave numbers for Fourier transform
    k = np.fft.fftfreq(n, d=dx) * 2 * np.pi  # Wave numbers in 1D
    laplacian_k = -k**2  # Fourier Laplacian operator

    for i in np.arange(ntimes): #Temporal loop
        rho_curr = rho[:, i]
        eta = np.sqrt(np.sum(free_energy_density(rho_curr) - C) * dx)  # Integral approximation via sum

        # 2. Compute mu
        f_rho_prime = F_der_contractive(rho_curr) - F_der_expansive(rho_curr)
        mu = -epsilon**2 * np.real(ifft(laplacian_k * fft(rho_curr))) + eta * f_rho_prime / eta
        
        # 3. Update rho in Fourier space
        mu_hat = fft(mu)
        rho_hat = fft(rho_curr)
        rho_hat_new = rho_hat - dt * M * laplacian_k * mu_hat
        rho_curr = np.real(ifft(rho_hat_new))
        
        # 4. Ensure physical bounds
        rho[:, i+1] = np.maximum(rho_curr, rho_min)

        if i % 10 == 0:
            X = (-1.0 + np.sqrt(1.0 + 2.0 * two_valence_delta * rho_curr)) / (two_valence_delta * rho_curr)
            grad = (np.roll(rho_curr, -1) - rho_curr) / dx
            F_bulk = rho_curr * (np.log(rho_curr) - 1.0) + B2 * rho_curr**2 + valence * rho_curr * (np.log(X) + 0.5 * (1.0 - X))
            F_interf = epsilon**2 / 2 * grad**2 # works in 1D
            F = np.sum(F_bulk + F_interf) * dx**3 / n

            print(f"{dt * i:.2g} {F}")
    
    with open("last_rho.dat", "w") as f:
        rho[:,-1].tofile(f, "\n")
        f.write("\n")
    

if __name__ == '__main__':
    main()