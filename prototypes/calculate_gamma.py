import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def estimate_surface_tension(rho_1, rho_2, M, B2_coeff, delta, kappa, sigma=1.0):
    """
    rho_1, rho_2: Coexisting vapor and liquid densities
    M: Valence (number of patches)
    epsilon: Bond energy (in units of kT)
    kappa: Influence parameter (gradient energy coeff)
    sigma: Particle diameter
    """
    
    # 1. Wertheim Association Helper
    def get_X(rho):
        # Solve quadratic: M*rho*delta*X^2 + X - 1 = 0
        return (-1.0 + np.sqrt(1.0 + 4.0 * M * delta * rho)) / (2.0 * M * delta * rho)

    def local_free_energy(rho):
        # Ideal + B_2 + Wertheim
        f_id = rho * (np.log(rho) - 1)
        f_B2 = B2_coeff * rho**2;
        
        X = get_X(rho)
        f_assoc = rho * M * (np.log(X) - X/2 + 0.5)
        return f_id + f_B2 + f_assoc

    # 2. Determine Coexistence Mu and P
    # (In a real scenario, these are calculated at the binodal)
    def get_mu_P(rho):
        X = get_X(rho)
        mu = np.log(rho) + 2.0 * B2_coeff * rho + M * np.log(X)
        P = mu * rho - local_free_energy(rho)
        return mu, P

    mu_1, P_1 = get_mu_P(rho_1)
    mu_2, P_2 = get_mu_P(rho_2)
    print(f"mu_1 = {mu_1}, P_1 = {P_1}")
    print(f"mu_2 = {mu_2}, P_2 = {P_2}")

    # 3. Define the ODE system for BVP
    # y[0] = rho, y[1] = d_rho/dx
    def ode_sys(x, y):
        rho = np.clip(y[0], 1e-10, 2.0)

        X = get_X(rho)
        mu = np.log(rho) + 2.0 * B2_coeff * rho + M * np.log(X)
        
        # d2_rho/dx2 = (mu_local - mu_coex) / kappa
        d2_rho = (mu - mu_1) / kappa
        return np.vstack((y[1], d2_rho))

    def bc(ya, yb):
        return np.array([ya[0] - rho_2, yb[0] - rho_1])

    # 4. Solve
    x = np.linspace(-10, 10, 100) # Grid in units of sigma
    y_guess = np.vstack((np.linspace(rho_2, rho_1, x.size), np.zeros(x.size)))
    
    res = solve_bvp(ode_sys, bc, x, y_guess) # res.y[0] is the density, res.y[1] is its gradient

    # 5. Integrate for Surface Tension: gamma = integral( kappa * (d_rho/dx)^2 )
    gamma = np.trapz(kappa * res.y[1]**2, res.x)
    
    return res.x, res.y[0], gamma

# Example Usage
x_coords, rho_profile, tension = estimate_surface_tension(
    rho_1=0.00191513, rho_2=1.521216, M=4, B2_coeff=0.5, delta=10, kappa=1
)

print(f"Estimated Surface Tension: {tension:.4f}")

plt.plot(x_coords, rho_profile)
plt.title("Interfacial Density Profile (Wertheim Theory)")
plt.xlabel("Position (x/σ)")
plt.ylabel("Density (ρ)")
plt.grid(True)
plt.show()
