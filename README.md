Monte carlo simulation of dark current in an irradiated semiconductor. The numerical constants for this simuatlion are:

# Physical Constants
q = 1.6e-19    # Elementary charge, C
k_B = 8.617e-5 # Boltzmann constant, eV/K
T0 = 300.0     # Reference temperature, K

# Material Properties (Silicon)
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
