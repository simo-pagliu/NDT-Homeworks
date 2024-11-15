import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
SECONDS_PER_DAY = 86400
DURATION_DAYS = 360
DURATION_SECONDS = DURATION_DAYS * SECONDS_PER_DAY
AVOGADRO = 6.022e23  # Avogadro's number (atoms/mol)

# Neutron flux parameters
peak_neutron_flux_fast = 6.1e15  # n/cm^2/s
peak_factors = np.array([0.572, 0.737, 0.868, 0.958, 0.983, 0.912, 0.802, 0.658, 0.498, 0.400])

# Neutron capture cross-sections (in cm^2)
sigma_f_Fe = 0.23e-24  # fast for Fe
sigma_f_Cr = 0.20e-24  # fast for Cr
sigma_f_Ni = 4.2e-24   # fast for Ni (assumed same for simplicity)
sigma_f_Ni58 = 4.4e-24 # fast for Ni-58
sigma_f_Ni59 = 1.3e-24 # fast for Ni-59

# Cladding composition (wt. %)
wt_Fe = 15.0
wt_Cr = 15.0
wt_Ni = 15.0
wt_Mo = 1.5
wt_Mn = 1.5
wt_Si = 0.9
wt_Ti = 0.4
wt_C = 0.09

# Material density for stainless steel (g/cm^3)
rho = 7.9  # g/cm³

# Atomic masses (g/mol)
M_Fe = 55.845
M_Cr = 52.00
M_Ni = 58.693
M_Mo = 95.95
M_Mn = 54.938
M_Si = 28.085
M_Ti = 47.867
M_C = 12.011

# Calculate atomic densities (atoms/cm³) for each element in cladding
def calc_atomic_density(wt_percent, molar_mass):
    """Calculate atomic density of an element in the cladding material."""
    return (wt_percent / molar_mass) * (rho * AVOGADRO)

n_Fe = calc_atomic_density(wt_Fe, M_Fe)
n_Cr = calc_atomic_density(wt_Cr, M_Cr)
n_Ni = calc_atomic_density(wt_Ni, M_Ni)
n_Mo = calc_atomic_density(wt_Mo, M_Mo)
n_Mn = calc_atomic_density(wt_Mn, M_Mn)
n_Si = calc_atomic_density(wt_Si, M_Si)
n_Ti = calc_atomic_density(wt_Ti, M_Ti)
n_C = calc_atomic_density(wt_C, M_C)

# Initial atomic densities
n_Ni58 = n_Ni * 0.68  # Ni-58 abundance
n_Ni59 = 0.0  # Ni-59 starts at zero

# Displacement cross-sections for DPA calculation (cm^2) 
sigma_dpa_Fe = 0.65e-24  # DPA cross-section for Fe
sigma_dpa_Cr = 0.7e-24   # DPA cross-section for Cr
sigma_dpa_Ni = 0.9e-24   # DPA cross-section for Ni

# Helium production rate adjustment (arbitrary scaling for clarity)
helium_scaling = 1e-10

def isotopic_changes(t, y):
    """Define the system of differential equations for isotope changes."""
    dNi58_dt = np.zeros(10)
    dNi59_dt = np.zeros(10)
    dFe_dt = np.zeros(10)
    dCr_dt = np.zeros(10)
    dHe_dt = np.zeros(10)  # helium concentration

    for i in range(10):
        # Neutron flux at each node
        phi_fast = peak_factors[i] * peak_neutron_flux_fast

        # Concentrations at node i
        Ni58 = y[i]
        Ni59 = y[i + 10]
        Fe = y[i + 20]
        Cr = y[i + 30]
        He = y[i + 40]  # helium concentration

        # Differential equations for isotopes
        dNi58_dt[i] = -sigma_f_Ni58 * phi_fast * Ni58
        dNi59_dt[i] = sigma_f_Ni58 * phi_fast * Ni58 - sigma_f_Ni59 * phi_fast * Ni59
        dFe_dt[i] = -sigma_f_Fe * phi_fast * Fe
        dCr_dt[i] = -sigma_f_Cr * phi_fast * Cr

        # Helium production from all sources
        dHe_dt[i] = (sigma_f_Fe * Fe + sigma_f_Cr * Cr + sigma_f_Ni * (Ni58 + Ni59)) * phi_fast * helium_scaling

    return np.concatenate((dNi58_dt, dNi59_dt, dFe_dt, dCr_dt, dHe_dt))

def helium_production_one_year(sigma, n, peak_factors, flux):
    """Calculate helium production for one year."""
    helium_prod = [sigma * n * factor * flux for factor in peak_factors]
    return sum(helium_prod) * DURATION_SECONDS  # total over 1 year

# Initial conditions vector
initial_conditions = np.concatenate((
    n_Ni58 * np.ones(10),
    n_Ni59 * np.ones(10),
    n_Fe * np.ones(10),
    n_Cr * np.ones(10),
    np.zeros(10)  # initial helium concentration
))

# Solve the differential equations
solution = solve_ivp(
    isotopic_changes,
    [0, DURATION_SECONDS],
    initial_conditions,
    method='RK45',
    dense_output=True
)

# Extract results
time_days = solution.t / SECONDS_PER_DAY
helium_concentration = solution.y[40:50, :]  # Extract helium concentration

# Calculate total helium production in atoms/cm^3 over 1 year
He_total_year_nodes = np.sum(helium_concentration[:, -1]) * DURATION_DAYS

# Calculate ppm of He for Fe (can be adjusted for Cr, Ni, etc.)
ppm_He = (He_total_year_nodes / n_Fe) * 1e6

# DPA calculations for each material
DPA_Fe = peak_neutron_flux_fast * sigma_dpa_Fe * DURATION_SECONDS / n_Fe
DPA_Cr = peak_neutron_flux_fast * sigma_dpa_Cr * DURATION_SECONDS / n_Cr
DPA_Ni = peak_neutron_flux_fast * sigma_dpa_Ni * DURATION_SECONDS / n_Ni

# Total DPA
total_DPA = DPA_Fe + DPA_Cr + DPA_Ni

# Ratio of ppmHe to DPA
ppm_He_to_DPA = ppm_He / total_DPA

# Output results
print(f"Total helium production in ppm: {ppm_He:.2f} ppm")
print(f"Total DPA: {total_DPA:.2e} DPA")
print(f"ppmHe / DPA ratio: {ppm_He_to_DPA:.2e}")

# Plot helium concentration over time for each node
plt.figure(figsize=(10, 6))
for i in range(10):
    plt.plot(time_days, helium_concentration[i], label=f'Node {i+1}')

plt.xlabel("Time (days)")
plt.ylabel("Helium Concentration (atoms/cm^3)")
plt.legend(loc="upper right")
plt.title("Helium Production Over Time at Each Node")
plt.show()

# Calculate helium production over 1 year for each element (atoms/cm^3)
He_Fe_year = helium_production_one_year(sigma_f_Fe, n_Fe, peak_factors, peak_neutron_flux_fast)
He_Cr_year = helium_production_one_year(sigma_f_Cr, n_Cr, peak_factors, peak_neutron_flux_fast)
He_Ni_year = helium_production_one_year(sigma_f_Ni, n_Ni, peak_factors, peak_neutron_flux_fast)
He_Ni58_year = helium_production_one_year(sigma_f_Ni58, n_Ni58, peak_factors, peak_neutron_flux_fast)
He_Ni59_year = helium_production_one_year(sigma_f_Ni59, n_Ni59, peak_factors, peak_neutron_flux_fast)

# Convert helium production to ppm (parts per million)
ppm_Fe_year = (He_Fe_year / n_Fe) * 1e6
ppm_Cr_year = (He_Cr_year / n_Cr) * 1e6
ppm_Ni_year = (He_Ni_year / n_Ni) * 1e6
ppm_Ni58_year = (He_Ni58_year / n_Ni58) * 1e6
ppm_Ni59_year = (He_Ni59_year / n_Ni59) * 1e6 if n_Ni59 > 0 else 0

# Display results for individual elements
print(f"Helium production due to Fe: {He_Fe_year:.3e} atoms/cm^3, {ppm_Fe_year:.3e} ppm")
print(f"Helium production due to Cr: {He_Cr_year:.3e} atoms/cm^3, {ppm_Cr_year:.3e} ppm")
print(f"Helium production due to Ni: {He_Ni_year:.3e} atoms/cm^3, {ppm_Ni_year:.3e} ppm")
print(f"Helium production due to Ni-58: {He_Ni58_year:.3e} atoms/cm^3, {ppm_Ni58_year:.3e} ppm")
print(f"Helium production due to Ni-59: {He_Ni59_year:.3e} atoms/cm^3, {ppm_Ni59_year:.3e} ppm")
