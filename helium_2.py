import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import ticker
from scipy.integrate import solve_ivp


# Constants
avogadro_number = 6.022e23  # atoms/mol
density_steel = 7.9  # g/cm³ (approximate for stainless steel)
#time_seconds = 365 * 24 * 3600  # 1 year in seconds
# Constants
time_seconds = 2 * 365 * 24 * 3600  # 2 years in seconds

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
    "B10_thermal": 3837e-27,  # 10B(n,α) thermal
    "B10_fast": 623e-27,  # 10B(n,α) fast
    "Fe_fast": 0.23e-27,  # Fe(n,α)
    "Cr_fast": 0.20e-27,  # Cr(n,α)
    "Ni58_fast": 4.2e-27,  # 58Ni(n,α)
    "Ni59_thermal": 13.4e-27,  # 59Ni(n,γ)
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
per_term_flux = 0
thermal_flux = per_term_flux* peak_flux  # Assumed fraction of thermal flux
node_numbers = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
peak_factors = np.array([0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.498])
flux_values = (1 - per_term_flux)* peak_flux * peak_factors  # Flux for each node

tot_at_He = 0
fuel_column_height = 0.850 #m
r_out_fuel = 0.00542 #m
volume_fuel_column = np.pi* (r_out_fuel/ 2) ** 2 * fuel_column_height # m³

# Initialize total helium production across all nodes
total_helium_production = 0
num_nodes = len(flux_values)  # Number of nodes
# Initialize arrays to store cumulative helium contributions at each time step
average_helium_contributions = np.zeros_like(t_eval)

for fast_flux in flux_values:
    # Bateman system of equations with helium contributions
    
    def Bateman_sys_with_He_contrib(t, state):
        N_Fe, N_Cr, N_58Ni, N_59Ni, N_extraNi, N_B10, N_He, He_Fe, He_Cr, He_Ni_fast, He_Ni_therm, He_B10 = state
        thermal_flux = 6.1e13  # HYPOTHESIS: Thermal flux value (given constant)
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

    # # Plot atomic concentrations over time
    # for i, label in enumerate(labels):
    #     plt.figure(figsize=(8, 5))
    #     plt.plot(sol_extended.t / (24 * 3600), sol_extended.y[i], color=colors[i], label=label)
    #     plt.title(f"Atomic Concentration Over Time: {label}")
    #     plt.ylabel("Concentration (atoms m⁻³)")
    #     plt.xlabel("Time (days)")
    #     plt.legend(loc="best")
    #     ax = plt.gca()
    #     ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2e"))
    #     plt.grid()
    #     plt.show()

    # Extract helium contributions
    He_contrib_Fe = sol_extended.y[7]
    He_contrib_Cr = sol_extended.y[8]
    He_contrib_Ni_fast = sol_extended.y[9]
    He_contrib_Ni_therm = sol_extended.y[10]
    He_contrib_B10 = sol_extended.y[11]

    # Calculate total helium contribution for this node
    helium_node_total = (
        He_contrib_Fe[-1]
        + He_contrib_Cr[-1]
        + He_contrib_Ni_fast[-1]
        + He_contrib_Ni_therm[-1]
        + He_contrib_B10[-1]
    )
    
    
    # Add to total helium production
    total_helium_production += helium_node_total
    # Calculate the average helium contributions over all nodes
    average_helium_contributions = total_helium_production / num_nodes

    
    # Plot the average helium contribution over time
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


    # Obtain the total He production 
    tot_at_He += He_contrib_Fe[-1] + He_contrib_Cr[-1] + He_contrib_Ni_fast[-1] + He_contrib_Ni_therm[-1] + He_contrib_B10[-1]


# Calculate the total initial number of atoms
total_initial_atoms = (
    N_Fe_0 
    + N_Cr_0 
    + N_Ni_0 
    + N_58Ni_0 
    + N_59Ni_0 
    + N_extraNi_0 
    + N_B10_0 
    + N_He_0
)


# Dimensions of the fuel column
fuel_pellet_outer_diameter = 5.42e-3  # meters
fuel_column_height = 0.85  # meters
# Calculate the volume of the fuel column
fuel_pin_volume = math.pi * 1e6 * (fuel_pellet_outer_diameter / 2) ** 2 * fuel_column_height
tot_at_He_V = tot_at_He / fuel_pin_volume
# Calculate average helium production per node
average_helium_per_node = tot_at_He_V / num_nodes
# Calculate average helium concentration in ppm
average_helium_concentration_ppm = (average_helium_per_node / total_initial_atoms) * 1e6
# Print the results

He_contrib_Fe_per_node = He_contrib_Fe[-1] / num_nodes
average_helium_concentration_Fe_ppm = (He_contrib_Fe_per_node/ total_initial_atoms) * 1e6
He_contrib_Cr_per_node = He_contrib_Cr[-1] / num_nodes
average_helium_concentration_Cr_ppm = (He_contrib_Cr_per_node / total_initial_atoms) * 1e6
He_contrib_Ni_therm_per_node = He_contrib_Ni_therm[-1] / num_nodes
average_helium_concentration_Ni_therm_ppm = (He_contrib_Ni_therm_per_node / total_initial_atoms) * 1e6
He_contrib_Ni_fast_per_node = He_contrib_Ni_fast[-1] / num_nodes
average_helium_concentration_Ni_fast_ppm = (He_contrib_Ni_fast_per_node / total_initial_atoms) * 1e6
He_contrib_B10_per_node = He_contrib_B10[-1] / num_nodes
average_helium_concentration_B10_ppm = (He_contrib_B10_per_node / total_initial_atoms) * 1e6

# Print the helium production contributions and final average concentration
print(f"He from Fe: {average_helium_concentration_Fe_ppm:.2f} ppm")
print(f"He from Cr: {average_helium_concentration_Cr_ppm:.2f} ppm")
print(f"He from Ni (fast): {average_helium_concentration_Ni_fast_ppm:.2f} ppm")
print(f"He from Ni (thermal): {average_helium_concentration_Ni_therm_ppm:.2f} ppm")
print(f"He from B-10: {average_helium_concentration_B10_ppm:.2f} ppm")
print(f"Avarage total helium concentration: {average_helium_concentration_ppm:.2f} ppm")
# Calculate average helium production
average_helium_production = tot_at_He / num_nodes
# Print the final average value
print(f"Total helium production: {tot_at_He:.2e} atoms")
print(f"Average helium production across all nodes: {average_helium_production:.2e} atoms")