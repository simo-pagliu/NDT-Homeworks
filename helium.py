import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.integrate import solve_ivp

# Constants
avogadro_number = 6.022e23  # atoms/mol
density_steel = 7.9  # g/cm³ (approximate for stainless steel)
time_seconds = 360 * 24 * 3600  # 1 year in seconds

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
    "B10_fast": 623e-27,  # 10B(n,α) fast
    "Fe_fast": 0.23e-27,  # Fe(n,α)
    "Cr_fast": 0.20e-27,  # Cr(n,α)
    "Ni58_fast": 4.2e-27,  # 58Ni(n,α)
}

# Isotopic abundance
isotopic_abundance = {"Ni58": 0.683, "B10": 0.198}

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
    [N_Fe_0, N_Cr_0, N_58Ni_0, N_59Ni_0, N_extraNi_0, N_B10_0, N_He_0, 0, 0, 0, 0]
)

# Axial nodalization data
node_numbers = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
peak_factors = np.array([0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.498])
peak_flux = 6.1e15  # n/cm²/s at the peak node
flux_values = peak_flux * peak_factors  # Flux for each node

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
    dN_59Ni_dt = 0 #cross_sections["Ni58_fast"] * flux * N_58Ni
    dN_extraNi_dt = 0 #-cross_sections["Ni58_fast"] * flux * N_extraNi
    dN_Fe_dt = -cross_sections["Fe_fast"] * flux * N_Fe
    dN_Cr_dt = -cross_sections["Cr_fast"] * flux * N_Cr
    dN_B10_dt = -cross_sections["B10_fast"] * flux * N_B10
    
    # Helium production contributions
    dHe_Fe = 0 #cross_sections["Fe_fast"] * flux * N_Fe
    dHe_Cr = 0 #cross_sections["Cr_fast"] * flux * N_Cr
    dHe_Ni_fast = 0 #cross_sections["Ni58_fast"] * flux * (N_58Ni + N_59Ni + N_extraNi)
    dHe_B10 = 0 #cross_sections["B10_fast"] * flux * N_B10
    dN_He_dt = dHe_Fe + dHe_Cr + dHe_Ni_fast + dHe_B10
    
    # Combine rates
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

# Solve the system for each node
results = []
for flux in [flux_values[1]]:
    sol = solve_ivp(
        lambda t, y: Bateman_sys_with_He_contrib(t, y, flux),
        t_span,
        N_0_extended,
        t_eval=t_eval,
        method="RK45",
    )
    results.append(sol)



# Plot helium concentration for each node
for node_index, sol in enumerate(results):
    sol_y = np.zeros_like(sol.y[0])
    for i in range(4):
        print(sol.y[i])
        sol_y =+ np.array(sol.y[i])

    plt.figure(figsize=(10, 6))
    plt.plot(
        sol.t / (24 * 3600),
        sol_y,
        label=f"Node {node_numbers[node_index]} (Flux={flux_values[node_index]:.2e})",
    )
    plt.title("Helium Concentration Over Time")
    plt.xlabel("Time (days)")
    plt.ylabel("Helium Concentration (atoms m⁻³)")
    plt.legend(loc="best")
    plt.grid()
    plt.show()
