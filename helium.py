import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
SECONDS_PER_DAY = 86400
DURATION_DAYS = 360
DURATION_SECONDS = DURATION_DAYS * SECONDS_PER_DAY

# Neutron flux parameters
peak_neutron_flux_fast = 6.1e15  # n/cm^2/s
peak_factors = np.array([0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.498])

# Neutron capture cross-sections (in cm^2)
sigma_th_Ni58 = 4.4e-24  # thermal for Ni-58
sigma_f_Ni58 = 4.2e-24   # fast for Ni-58 (assumed same for simplicity)
sigma_f_Ni59 = 1.3e-24   # fast for Ni-59
sigma_f_Fe = 0.23e-24    # fast for Fe
sigma_f_Cr = 0.20e-24    # fast for Cr

# Initial concentrations (arbitrary units)
Ni58_initial = 1.0
Ni59_initial = 0.0
Fe_initial = 1.0
Cr_initial = 1.0

# Define the system of differential equations
def isotopic_changes(t, y):
    dNi58_dt = np.zeros(10)
    dNi59_dt = np.zeros(10)
    dFe_dt = np.zeros(10)
    dCr_dt = np.zeros(10)

    for i in range(10):
        # Neutron flux at each node
        phi_fast = peak_factors[i] * peak_neutron_flux_fast

        # Concentrations at node i
        Ni58 = y[i]
        Ni59 = y[i + 10]
        Fe = y[i + 20]
        Cr = y[i + 30]

        # Differential equations
        dNi58_dt[i] = -sigma_th_Ni58 * phi_fast * Ni58 - sigma_f_Ni58 * phi_fast * Ni58
        dNi59_dt[i] = sigma_th_Ni58 * phi_fast * Ni58 - sigma_f_Ni59 * phi_fast * Ni59
        dFe_dt[i] = -sigma_f_Fe * phi_fast * Fe
        dCr_dt[i] = -sigma_f_Cr * phi_fast * Cr

    return np.concatenate((dNi58_dt, dNi59_dt, dFe_dt, dCr_dt))

# Initial conditions vector
initial_conditions = np.concatenate((
    Ni58_initial * np.ones(10),
    Ni59_initial * np.ones(10),
    Fe_initial * np.ones(10),
    Cr_initial * np.ones(10)
))

# Solve the differential equations
solution = solve_ivp(
    isotopic_changes,
    [0, DURATION_SECONDS],
    initial_conditions,
    method='RK45',
    dense_output=True
)

# Plotting the results
time_days = solution.t / SECONDS_PER_DAY
for i in range(10):
    plt.plot(time_days, solution.y[i], label=f'Ni-58 Node {i+1}')
    plt.plot(time_days, solution.y[i+10], label=f'Ni-59 Node {i+1}')
    plt.plot(time_days, solution.y[i+20], label=f'Fe Node {i+1}')
    plt.plot(time_days, solution.y[i+30], label=f'Cr Node {i+1}')

plt.xlabel("Time (days)")
plt.ylabel("Concentration")
plt.legend(loc="upper right", bbox_to_anchor=(1.2, 1))
plt.title("Isotopic Concentrations Over Time at Each Node")
plt.show()

