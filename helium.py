import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.integrate import solve_ivp

# Constants
avogadro_number = 6.022e23  # atoms/mol
density_steel = 7.9  # g/cm³ (approximate for stainless steel)
time_seconds = 365 * 24 * 3600  # 1 year in seconds

# Cladding composition (wt.%)
composition = {
    "Fe": 100 - (15.0 + 15.0 + 1.5 + 1.5 + 0.9 + 0.4 + 0.09),  # Balance is Iron
    "Cr": 15.0,
    "Ni": 15.0,
    "B": 0.006,  # Boron content in wt.% (ppm equivalent: 60 ppm)
}
molar_masses = {
    "Fe": 55.847,
    "Cr": 51.996,
    "Ni": 58.690,
    "B": 10.811,
}

# Cross-sections (cm²)
cross_sections = {
    "B10_fast": 623e-24,  # 10B(n,α) fast
    "Fe_fast": 0.23e-24,  # Fe(n,α)
    "Cr_fast": 0.20e-24,  # Cr(n,α)
    "Ni58_fast": 4.2e-24,  # 58Ni(n,α)
}

# Isotopic abundance
isotopic_abundance = {"Ni58": 0.683, "B10": 0.198}

# Neutron flux data for each node
node_centers = np.array([42.5, 127.5, 212.5, 297.5, 382.5, 467.5, 552.5, 637.5, 722.5, 807.5])  # mm
peak_factors = np.array([0.572, 0.737, 0.868, 0.958, 1.0, 0.983, 0.912, 0.802, 0.658, 0.498])
peak_flux = 6.1e15  # n/cm²/s at peak node
Flux_values = peak_flux * peak_factors  # Flux at each node

# Height differences for volumes
node_heights = np.diff(node_centers, prepend=0)  # Add a 0 at the start for full range
volumes = node_heights  # Volume proportional to height difference (assuming constant cross-sectional area)

# Helper function for initial concentrations
def concentration(wt_percent, density, molar_mass):
    return (wt_percent * density * avogadro_number) / molar_mass

# Initial concentrations
N_Fe_0 = concentration(composition["Fe"], density_steel, molar_masses["Fe"])
N_Cr_0 = concentration(composition["Cr"], density_steel, molar_masses["Cr"])
N_Ni_0 = concentration(composition["Ni"], density_steel, molar_masses["Ni"])
N_58Ni_0 = isotopic_abundance["Ni58"] * N_Ni_0
N_59Ni_0 = 0  # Initial 59Ni is 0
N_extraNi_0 = N_Ni_0 * (1 - isotopic_abundance["Ni58"])
N_B10_0 = concentration(composition["B"], density_steel, molar_masses["B"]) * isotopic_abundance["B10"]
N_He_0 = 0  # Initial helium concentration

# Extended state vector
N_0_extended = np.array(
    [N_Fe_0, N_Cr_0, N_58Ni_0, N_59Ni_0, N_extraNi_0, N_B10_0, N_He_0, 0, 0, 0, 0, 0]
)

# Time parameters
t_0 = 0
t_f = time_seconds
t_span = (t_0, t_f)
t_eval = np.linspace(t_0, t_f, 1000)

# Bateman system of equations with helium contributions
def Bateman_sys_with_He_contrib(t, state, flux):
    N_Fe, N_Cr, N_58Ni, N_59Ni, N_extraNi, N_B10, N_He, He_Fe, He_Cr, He_Ni_fast, He_B10 = state

    # Reaction rates
    dN_58Ni_dt = -cross_sections["Ni58_fast"] * flux * N_58Ni
    dN_59Ni_dt = cross_sections["Ni58_fast"] * flux * N_58Ni
    dN_extraNi_dt = -cross_sections["Ni58_fast"] * flux * N_extraNi
    dN_Fe_dt = -cross_sections["Fe_fast"] * flux * N_Fe
    dN_Cr_dt = -cross_sections["Cr_fast"] * flux * N_Cr
    dN_B10_dt = -cross_sections["B10_fast"] * flux * N_B10

    # Helium production contributions
    dHe_Fe = cross_sections["Fe_fast"] * flux * N_Fe
    dHe_Cr = cross_sections["Cr_fast"] * flux * N_Cr
    dHe_Ni_fast = cross_sections["Ni58_fast"] * flux * (N_58Ni + N_59Ni + N_extraNi)
    dHe_B10 = cross_sections["B10_fast"] * flux * N_B10
    dN_He_dt = dHe_Fe + dHe_Cr + dHe_Ni_fast + dHe_B10

    return [
        dN_Fe_dt,
        dN_Cr_dt,
        dN_58Ni_dt,
        dN_59Ni_dt,
        dN_extraNi_dt,
        dN_B10_dt,
        dN_He_dt,
        dHe_Fe,
        dHe_Cr,
        dHe_Ni_fast,
        dHe_B10,
    ]

# Solve the system for each node and store results
helium_concentrations = []
for flux in Flux_values:
    sol = solve_ivp(
        Bateman_sys_with_He_contrib,
        t_span,
        N_0_extended,
        t_eval=t_eval,
        method="RK45",
        args=(flux,),
    )
    helium_concentrations.append(sol.y[6])  # N_He is at index 6

# Calculate weighted average helium concentration
helium_concentrations = np.array(helium_concentrations)
weighted_average_helium = np.sum(helium_concentrations[:, -1] * volumes) / np.sum(volumes)

print(f"Weighted Average Helium Concentration: {weighted_average_helium:.2e} atoms/m³")