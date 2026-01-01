# SIRR Model Simulation (with Vital Dynamics)

This repository contains a Python implementation of an extended epidemiological model: **SIRR (Susceptible - Infected - Recovered - Reinfected)**. Unlike standard SIR models, this simulation accounts for reinfection and vital dynamics (birth and death rates), providing a more complex view of disease propagation over time.

## Project Overview

The simulation models the flow of a population through four distinct states:
1.  **Susceptible (S):** Individuals who can catch the disease.
2.  **Infected (I):** Individuals currently infected.
3.  **Recovered (R):** Individuals who have recovered from the disease.
4.  **Reinfected (Re):** Individuals who have lost immunity and become reinfected (or are tracked as a specific reinfection compartment).

### Key Features
* **Differential Equation Modeling:** Uses `scipy.integrate.odeint` to solve the system of Ordinary Differential Equations (ODEs).
* **Vital Dynamics:** Includes parameters for birth rates and death rates, allowing for population fluctuation over long periods.
* **Parameter Estimation:** Generates noisy synthetic data and uses `scipy.optimize.curve_fit` to estimate the original epidemiological parameters ($\beta$, $\gamma$, $\delta$, etc.).
* **Visualization:** Plots the population dynamics of all compartments over time using `matplotlib`.

## Mathematical Model

The system is defined by the following set of differential equations:

$$\frac{dS}{dt} = -\beta SI + \mu_{birth} (S+I+R+R_e) - \mu_{death} S$$
$$\frac{dI}{dt} = \beta SI - \gamma I + \delta R - \mu_{death} I$$
$$\frac{dR}{dt} = \gamma I - \delta R - \mu_{death} R$$
$$\frac{dR_e}{dt} = \delta R - \mu_{death} R_e$$

*Where:*
* $\beta$: Infection rate
* $\gamma$: Recovery rate
* $\delta$: Reinfection/Immunity loss rate
* $\mu_{birth/death}$: Vital dynamic rates

## Installation & Usage

### Prerequisites
Ensure you have Python installed along with the following libraries:

```bash
pip install numpy scipy matplotlib
