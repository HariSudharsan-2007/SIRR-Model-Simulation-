import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define the SIRR Model
def sirr_model(y, t, beta, gamma, delta, birthrate, deathrate):
    S, I, R, Re = y
    
    # Population flow equations
    dSdt = -beta * S * I + birthrate * (S + I + R + Re) - deathrate * S
    dIdt = beta * S * I - gamma * I + delta * R - deathrate * I
    
    dRdt = gamma * I - delta * R - deathrate * R
    dRedt = delta * R - deathrate * Re
    
    return [dSdt, dIdt, dRdt, dRedt]

# Helper Functions for Solving
def solve_sirr_full(t, beta, gamma, delta, birthrate, deathrate, S0, I0, R0, Re0):
    y0 = [S0, I0, R0, Re0]
    ret = odeint(sirr_model, y0, t, args=(beta, gamma, delta, birthrate, deathrate))
    return ret

def solve_sirr_infected_only(t, beta, gamma, delta, birthrate, deathrate, S0, I0, R0, Re0):
    # Wrapper function to return only the Infected (I) count
    # Used specifically for curve_fit optimization
    return solve_sirr_full(t, beta, gamma, delta, birthrate, deathrate, S0, I0, R0, Re0)[:, 1]

# Simulation and Parameter Estimation

# Generate Simulated Data
t_data = np.linspace(0, 50, 365)
true_params = {
    'beta': 0.3,
    'gamma': 0.05,
    'delta': 0.02,
    'birthrate': 0.01,
    'deathrate': 0.01
}

# Initial conditions for data generation
S0_data, I0_data, R0_data, Re0_data = 1000, 10, 0, 0

I_data = solve_sirr_infected_only(t_data, **true_params, S0=S0_data, I0=I0_data, R0=R0_data, Re0=Re0_data)

# Add Gaussian noise to simulate measurement error
noise = np.random.normal(0, 5, size=len(I_data))
I_data_noisy = I_data + noise

# Fit the Model to the Noisy Data
popt, _ = curve_fit(
    lambda t, beta, gamma, delta, birthrate, deathrate: solve_sirr_infected_only(
        t, beta, gamma, delta, birthrate, deathrate, S0=S0_data, I0=I0_data, R0=R0_data, Re0=Re0_data
    ),
    t_data,
    I_data_noisy,
    bounds=(0, [1.0, 1.0, 1.0, 0.1, 0.1])
)

beta_est, gamma_est, delta_est, birthrate_est, deathrate_est = popt
print(f"Estimated Parameters:")
print(f"Beta: {beta_est:.4f}")
print(f"Gamma: {gamma_est:.4f}")
print(f"Delta: {delta_est:.4f}")
print(f"Birth Rate: {birthrate_est:.4f}")
print(f"Death Rate: {deathrate_est:.4f}")

# Visualization
# Generate smooth curves using the ESTIMATED parameters
t_smooth = np.linspace(0, 50, 200)
SIRR_curves = solve_sirr_full(t_smooth, beta_est, gamma_est, delta_est, birthrate_est, deathrate_est, S0_data, I0_data, R0_data, Re0_data)

S_curve = SIRR_curves[:, 0]
I_curve = SIRR_curves[:, 1]
R_curve = SIRR_curves[:, 2]
Re_curve = SIRR_curves[:, 3]

plt.figure(figsize=(10, 6))
plt.plot(t_smooth, S_curve, label='Susceptible (S)', color='blue')
plt.plot(t_smooth, I_curve, label='Infected (I)', color='red')
plt.plot(t_smooth, R_curve, label='Recovered (R)', color='green')
plt.plot(t_smooth, Re_curve, label='Reinfected (Re)', color='purple')

plt.title("SIRR Model Simulation (Including Birthrate and Deathrate)")
plt.xlabel("Time")
plt.ylabel("Population")
plt.legend()
plt.grid(True)
plt.show()