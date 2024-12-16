import numpy as np
import matplotlib.pyplot as plt
import random

# Physical Constants
q = 1.6e-19    # Elementary charge, C
k_B = 8.617e-5 # Boltzmann constant, eV/K
T0 = 300.0     # Reference temperature, K

# Material Properties (Silicon example)
E_g0 = 1.12     # Bandgap energy at T0=300K, eV
N_c0 = 2.8e19   # Eff. density of states conduction band at 300K, cm^-3
N_v0 = 1.04e19  # Eff. density of states valence band at 300K, cm^-3
alpha = 4.73e-4 # Varshni parameter (eV/K)
beta = 636      # Varshni parameter (K)

# Defect Properties (SRH)
E_DL = 0.6          # Deep level energy above valence band, eV
sigma_e = 1e-15     # Capture cross-section for electrons, cm^2
sigma_h = 1e-15     # Capture cross-section for holes, cm^2
v_e = 1e7           # Thermal velocity of electrons, cm/s
v_h = 1e7           # Thermal velocity of holes, cm/s
N_DL_initial = 1e15 # Initial defect density, cm^-3

# Auger Recombination Parameters
C_Auger = 1e-31     # Auger recombination coefficient, cm^6/s

# Geometry and Timescale
W = 1e-4   # Depletion width, cm
A = 1.0    # Area, cm^2
Volume = W * A  # Volume, cm^3

# Mitigation and Time Parameters
f_anneal = 1.0   # Annealing factor (no annealing if =1)
n_splits = 4     # For splitting the device into sub-sections
tau_anneal = 200  # Annealing time scale (s)
G = 1e10          # Generation rate, cm^-3 s^-1

# Dose parameters
dose_scale = 1e11 # Scale for dose exponential

# Equations

def N_c(T):
    #Compute effective density of states in conduction band at temperature T: (T/T0)^(3/2).

    return N_c0 * (T/T0)**1.5

def N_v(T):

    #Compute effective density of states in valence band at temperature T:  as (T/T0)^(3/2).

    return N_v0 * (T/T0)**1.5

def E_g(T):
    #Compute the bandgap E_g at temperature T using the Varshni equation: E_g(T) = E_g0 - (alpha * T^2 / (T + beta))

    return E_g0 - (alpha*T**2/(T+beta))

def n_i(T):
    #Compute intrinsic carrier concentration n_i at temperature T: n_i(T) = sqrt(N_c(T)*N_v(T)) * exp(-E_g(T)/(2*k_B*T))

    Eg = E_g(T)
    return np.sqrt(N_c(T)*N_v(T))*np.exp(-Eg/(2*k_B*T))

def apply_splitting_and_annealing(N_DL, n_split, factor):
    #Divides by n_split and multipleis by factor. Returns N_DL_eff.

    return N_DL/(n_split)*factor

def calculate_SRH_rate(n, p, N_DL_eff, E_DL, v_e, sigma_e, v_h, sigma_h, T):
    #Calculate the Shockley-Read-Hall (SRH) recombination rate: U_SRH = (v_e*sigma_e*N_DL_eff*(n*p - n_i^2)) / (n + p + 2*n_i*cosh((E_DL - E_MB)/(k_B*T)))

    E_MB = E_g(T)/2
    ni_T = n_i(T)
    term1 = n*p - ni_T**2
    denominator = n + p + 2*ni_T*np.cosh((E_DL - E_MB)/(k_B*T))
    U = (v_e * sigma_e * N_DL_eff * term1)/denominator
    return U

def calculate_Auger_rate(n, p, C_Auger):
    #Calculate Auger recombination rate: R_Auger = C_Auger * n * p * (n + p)

    return C_Auger*n*p*(n+p)

def total_recombination_rate(n, p, N_DL_eff, E_DL, v_e, sigma_e, v_h, sigma_h, T, C_Auger):
    #Compute total recombination rate as sum of SRH and Auger: R_total = R_SRH + R_Auger

    R_SRH = calculate_SRH_rate(n, p, N_DL_eff, E_DL, v_e, sigma_e, v_h, sigma_h, T)
    R_Auger = calculate_Auger_rate(n, p, C_Auger)
    return R_SRH + R_Auger

def generate_carriers(n, p, G, dt):
    #Generate carriers over time step dt: dn = dp = G * dt

    dn = G*dt
    dp = G*dt
    return n+dn, p+dp

def monte_carlo_iteration(n, p, T, N_DL_eff, G, Delta_t=1e-7):
    """
    Perform one simulation iteration:
    - Generate carriers
    - total recombination rate
    - Recombine carriers
    - Compute dark current J_d = q * R_total * W
    """
    n, p = generate_carriers(n, p, G, Delta_t)

    R_total = total_recombination_rate(n, p, N_DL_eff, E_DL, v_e, sigma_e, v_h, sigma_h, T, C_Auger)

    dR = R_total * Delta_t
    n = max(n - dR, 0)
    p = max(p - dR, 0)

    J_d = q * R_total * W
    return n, p, J_d


# Simulation Conditions
num_iterations = 100
Delta_t_sim = 1e-7 # stable time step

def simulate_vs_time(T=300.0):
    """
    Simulate dark current over time (0 to 1000 s) under annealing. N_DL_eff(t) = (N_DL_initial/n_splits)*exp(-t/tau_anneal).
    Returns array of (t, J_d).
    """

    times = np.linspace(0,1000,50)
    n = p = n_i(T)
    results_t = []

    for t in times:
        factor_anneal = np.exp(-t/tau_anneal)
        N_DL_eff = (N_DL_initial/n_splits)*factor_anneal
        n, p, J_d = monte_carlo_iteration(n, p, T, N_DL_eff, G, Delta_t=Delta_t_sim)
        results_t.append((t,J_d))
    return np.array(results_t)

def simulate_vs_dose(T=300.0):
    """
    Simulate dark current as a function of dose.
    N_DL(dose) = (N_DL_initial * exp(dose/dose_scale))/n_splits.
    Returns array of (dose, J_d).
    """

    doses = np.linspace(0,5e11,50)
    n = p = n_i(T)
    results_d = []

    for dose in doses:
        factor_dose = np.exp(dose/dose_scale)
        N_DL_eff = (N_DL_initial*factor_dose)/n_splits
        n, p, J_d = monte_carlo_iteration(n, p, T, N_DL_eff, G, Delta_t=Delta_t_sim)
        results_d.append((dose,J_d))
    return np.array(results_d)

# Run simulations and plot

#Dark Current vs Temperature
N_DL_eff_base = N_DL_initial
split_list = [1, 4, 16]  # no split, 2x2 split, and 4x4 split
temperatures_C = np.arange(0, 200, 20)
temperatures_K = temperatures_C + 273.15

fig, ax = plt.subplots(figsize=(8,6))

for splits, color in zip(split_list, ['blue','green','red']):
    J_data = []
    for Temp in temperatures_K:
        # Apply splitting (no annealing)
        N_DL_eff = apply_splitting_and_annealing(N_DL_eff_base, splits, 1.0)
        n = p = n_i(Temp)
        J_vals = []
        for _ in range(num_iterations):
            n, p, J_d = monte_carlo_iteration(n, p, Temp, N_DL_eff, G, Delta_t=Delta_t_sim)
            J_vals.append(J_d)
        J_avg = np.mean(J_vals[int(num_iterations/2):])
        J_data.append(J_avg)

    ax.plot(temperatures_C, J_data, marker='o', color=color, label=f'{splits} splits')

ax.set_yscale('log')
ax.set_xlabel('Temperature (°C)')
ax.set_ylabel('Dark Current (A/cm²)')
ax.set_title('Dark Current vs Temperature for Different Splits (with Auger)')
ax.grid(True, which='both', ls='--')
ax.legend()
plt.tight_layout()
plt.show()

#Dark Current vs Time
results_t = simulate_vs_time(T=300.0) 
plt.figure()
plt.plot(results_t[:,0], results_t[:,1], marker='s', label='DC vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Dark Current (A/cm²)')
plt.yscale('log')
plt.grid(True, which='both', ls='--')
plt.title('Dark Current vs Time (Annealing)')
plt.legend()
plt.tight_layout()
plt.show()

#Dark Current vs Dose
results_d = simulate_vs_dose(T=300.0)
plt.figure()
plt.plot(results_d[:,0], results_d[:,1], marker='^', label='DC vs Dose')
plt.xlabel('Dose (arbitrary units)')
plt.ylabel('Dark Current (A/cm²)')
plt.yscale('log')
plt.grid(True, which='both', ls='--')
plt.title('Dark Current vs Received Dose')
plt.legend()
plt.tight_layout()
plt.show()