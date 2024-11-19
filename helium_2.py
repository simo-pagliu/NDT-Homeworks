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
    "B10_thermal": 3837e-24,  # 10B(n,α) thermal
    "B10_fast": 623e-24,  # 10B(n,α) fast
    "Fe_fast": 0.23e-24,  # Fe(n,α)
    "Cr_fast": 0.20e-24,  # Cr(n,α)
    "Ni58_fast": 4.2e-24,  # 58Ni(n,α)
    "Ni59_thermal": 13.4e-24,  # 59Ni(n,γ)
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
    [N_Fe_0, N_Cr_0, N_58Ni_0, N_59Ni_0, N_extraNi_0, N_B10_0, N_He_0, 0, 0, 0, 0, 0]
)

# Time parameters
t_0 = 0
t_f = time_seconds
t_span = (t_0, t_f)
t_eval = np.linspace(t_0, t_f, 1000)

# Neutron flux data
peak_flux = 6.1e15  # n/cm²/s at the peak node
thermal_flux = 0.8 * peak_flux  # Assumed fraction of thermal flux
fast_flux = 0.2 * peak_flux  # Assumed fraction of fast flux

# Bateman system of equations with helium contributions
def Bateman_sys_with_He_contrib(t, state):
    N_Fe, N_Cr, N_58Ni, N_59Ni, N_extraNi, N_B10, N_He, He_Fe, He_Cr, He_Ni_fast, He_Ni_therm, He_B10 = state
    
    # Reaction rates
    dN_58Ni_dt = -cross_sections["Ni58_fast"] * fast_flux * N_58Ni
    dN_59Ni_dt = (
        cross_sections["Ni58_fast"] * fast_flux * N_58Ni
        - cross_sections["Ni59_thermal"] * thermal_flux * N_59Ni
    )
    dN_extraNi_dt = -cross_sections["Ni58_fast"] * fast_flux * N_extraNi
    dN_Fe_dt = -cross_sections["Fe_fast"] * fast_flux * N_Fe
    dN_Cr_dt = -cross_sections["Cr_fast"] * fast_flux * N_Cr
    dN_B10_dt = -(
        cross_sections["B10_thermal"] * thermal_flux
        + cross_sections["B10_fast"] * fast_flux
    ) * N_B10
    
    # Helium production contributions
    dHe_Fe = cross_sections["Fe_fast"] * fast_flux * N_Fe
    dHe_Cr = cross_sections["Cr_fast"] * fast_flux * N_Cr
    dHe_Ni_fast = cross_sections["Ni58_fast"] * fast_flux * (N_58Ni + N_59Ni + N_extraNi)
    dHe_Ni_therm = cross_sections["Ni59_thermal"] * thermal_flux * N_59Ni
    dHe_B10 = (
        cross_sections["B10_thermal"] * thermal_flux * N_B10
        + cross_sections["B10_fast"] * fast_flux * N_B10
    )
    dN_He_dt = dHe_Fe + dHe_Cr + dHe_Ni_fast + dHe_Ni_therm + dHe_B10
    
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
        dHe_Ni_therm,
        dHe_B10,
    ]

# Solve the system
sol_extended = solve_ivp(
    Bateman_sys_with_He_contrib,
    t_span,
    N_0_extended,
    t_eval=t_eval,
    method="RK45",
)

# Labels and colors for the plots
labels = ["[Fe]", "[Cr]", "[Ni-58]", "[Ni-59]", "[Ni-extra]", "[B-10]", "[He]"]
colors = ["c", "b", "g", "y", "k", "m", "r"]

# Plot atomic concentrations over time
for i, label in enumerate(labels):
    plt.figure(figsize=(8, 5))
    plt.plot(sol_extended.t / (24 * 3600), sol_extended.y[i], color=colors[i], label=label)
    plt.title(f"Atomic Concentration Over Time: {label}")
    plt.ylabel("Concentration (atoms m⁻³)")
    plt.xlabel("Time (days)")
    plt.legend(loc="best")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2e"))
    plt.grid()
    plt.show()

# Extract helium contributions
He_contrib_Fe = sol_extended.y[7]
He_contrib_Cr = sol_extended.y[8]
He_contrib_Ni_fast = sol_extended.y[9]
He_contrib_Ni_therm = sol_extended.y[10]
He_contrib_B10 = sol_extended.y[11]

# Plot helium contributions over time
plt.figure(figsize=(10, 6))
plt.plot(sol_extended.t / (24 * 3600), He_contrib_Fe, label="He from Fe", color="c")
plt.plot(sol_extended.t / (24 * 3600), He_contrib_Cr, label="He from Cr", color="b")
plt.plot(sol_extended.t / (24 * 3600), He_contrib_Ni_fast, label="He from Ni (fast)", color="g")
plt.plot(sol_extended.t / (24 * 3600), He_contrib_Ni_therm, label="He from Ni (thermal)", color="y")
plt.plot(sol_extended.t / (24 * 3600), He_contrib_B10, label="He from B-10", color="m")
plt.title("Helium Contribution from Different Elements Over Time")
plt.xlabel("Time (days)")
plt.ylabel("Helium Production Contribution (atoms)")
plt.legend(loc="best")
plt.grid()
plt.show()
