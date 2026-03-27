import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def estimate_surface_tension(rho_1, rho_2, M, B2_coeff, delta, kappa):
    """
    rho_1, rho_2: Coexisting vapor and liquid densities
    M: Valence (number of patches)
    epsilon: Bond energy (in units of kT)
    kappa: Influence parameter (gradient energy coeff)
    """
    
    # 1. Wertheim Association Helper
    def get_X(rho):
        # Stable version of the quadratic solution
        # Avoids division by rho
        return 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * M * delta * rho))

    def beta_f(rho):
        # Ideal + B_2 + Wertheim
        f_id = rho * (np.log(rho) - 1)
        f_B2 = B2_coeff * rho**2;
        
        X = get_X(rho)
        f_assoc = rho * M * (np.log(X) - X/2 + 0.5)
        return f_id + f_B2 + f_assoc
    
    def get_mu(rho):
        X = get_X(rho)
        mu = np.log(rho) + 2.0 * B2_coeff * rho + M * np.log(X)
        return mu

    # 2. Determine Coexistence Mu and P
    # (In a real scenario, these are calculated at the binodal)
    def get_mu_P(rho):
        mu = get_mu(rho)
        P = mu * rho - beta_f(rho)
        return mu, P

    mu_1, P_1 = get_mu_P(rho_1)
    mu_2, P_2 = get_mu_P(rho_2)
    print(f"mu_1 = {mu_1}, P_1 = {P_1}")
    print(f"mu_2 = {mu_2}, P_2 = {P_2}")

    # 3. Define the ODE system for BVP
    # y[0] = rho, y[1] = d_rho/dx
    def ode_sys(x, y):
        rho = np.clip(y[0], 1e-10, 2.0)

        mu = get_mu(rho)
        
        # d2_rho/dx2 = (mu_local - mu_coex) / kappa
        d2_rho = (mu - mu_1) / (2.0 * kappa)
        return np.vstack((y[1], d2_rho))

    def bc(ya, yb):
        return np.array([ya[0] - rho_2, yb[0] - rho_1])

    # 4. Solve
    x = np.linspace(-50, 50, 1000) # Grid in units of sigma
    #y_guess = np.vstack((np.linspace(rho_2, rho_1, x.size), np.zeros(x.size)))

    # this value should be such that rho_guess is equal to rho_1 and rho_2 at the boundaries
    width_guess = 1.0
    rho_guess = (rho_2 + rho_1)/2 - (rho_2 - rho_1)/2 * np.tanh(x / width_guess)
    grad_guess = -(rho_2 - rho_1)/(2 * width_guess) * (1 / np.cosh(x / width_guess)**2)
    y_guess = np.vstack((rho_guess, grad_guess))
    
    # res.y[0] is the density, res.y[1] is its gradient
    res = solve_bvp(ode_sys, bc, x, y_guess, max_nodes=5000) # max_nodes should be chosen carefully

    # 5. Integrate for Surface Tension: gamma = integral( kappa * (d_rho/dx)^2 )
    gamma = np.trapz(2.0 * kappa * res.y[1]**2, res.x)
    
    return res.x, res.y[0], gamma

# Example Usage
x_coords, rho_profile, tension = estimate_surface_tension(
    rho_2=0.00191513, rho_1=1.521216, M=4, B2_coeff=0.5, delta=10, kappa=1
)

print(f"Estimated Surface Tension: {tension:.4f}")

np.savetxt("equilibrated_profile.dat", np.c_[x_coords, rho_profile])

plt.plot(x_coords, rho_profile)
plt.title("Interfacial Density Profile (Wertheim Theory)")
plt.xlabel("Position (x/σ)")
plt.ylabel("Density (ρ)")
plt.grid(True)
plt.show()
