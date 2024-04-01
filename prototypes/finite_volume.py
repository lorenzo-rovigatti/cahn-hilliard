# some of the code has been taken from https://github.com/sergiopperez/Cahn_Hilliard_Finite_Volume

import numpy as np
from scipy import optimize

def F_der_contractive(rho):
    B2 = 10600
    F_der = np.log(rho) + 2.0 * B2 * rho

    return F_der

def F_der_expansive(rho):
    M = 4
    Delta = 59769.7978988913
    two_valence_delta = 2 * M * Delta
    X = (-1.0 + np.sqrt(1.0 + 2.0 * two_valence_delta * rho)) / (two_valence_delta * rho)

    F_der = -M * np.log(X)

    return F_der

def to_solve(rho_curr, rho_next, epsilon, dx, dt, M, n):
    F_der_con = F_der_contractive(rho_next)
    F_der_exp = F_der_expansive(rho_curr)

    Lap = np.zeros(np.shape(rho_curr)[0])
    Lap[1:-1] = epsilon**2*(rho_curr[0:-2]-2*rho_curr[1:-1]+rho_curr[2:]+rho_next[0:-2]-2*rho_next[1:-1]+rho_next[2:])/dx**2/2.
    Lap[0] = epsilon**2*(rho_curr[-1]-2*rho_curr[0]+rho_curr[1]+rho_next[-1]-2*rho_next[0]+rho_next[1])/dx**2/2.
    Lap[-1] = epsilon**2*(rho_curr[-2]-2*rho_curr[-1]+rho_curr[0]+rho_next[-2]-2*rho_next[-1]+rho_next[0])/dx**2/2.

    uhalf = -(np.roll(F_der_con - F_der_exp - Lap, 1) - F_der_con + F_der_exp + Lap) / dx

    positive_indices_u = uhalf > 0
    negative_indices_u = uhalf < 0

    uhalfplus = np.zeros(n)
    uhalfminus = np.zeros(n)
    uhalfplus[positive_indices_u] = uhalf[positive_indices_u]
    uhalfminus[negative_indices_u] = uhalf[negative_indices_u]

    # Compute (n+1) fluxes, including no-flux boundary conditions
    Fhalf = np.zeros(n + 1)

    # 1st order
    Fhalf = uhalfplus * M + uhalfminus * M

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
    t = np.zeros(ntimes+1) # Time vector

    for i in np.arange(ntimes): #Temporal loop
        
        rho[:, i+1], infodict, ier, mesg = optimize.fsolve(lambda rhon: to_solve(rho[:,i], rhon, epsilon, dx, dt, M, n), rho[:,i], full_output = True)
        
        t[i+1] = t[i]+dt
        
        print('--------------------')
        print('Time: ', t[i])
        print(['L1 norm of the difference between the new and old state: ', np.linalg.norm(rho[:,i+1]-rho[:,i],1)])
    
    
    with open("last_rho.dat", "w") as f:
        rho[:,-1].tofile(f, "\n")
        f.write("\n")
 
    # Compute free energy in time
    F = 0
    

if __name__ == '__main__':
    main()