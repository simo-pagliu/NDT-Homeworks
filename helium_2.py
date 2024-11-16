import numpy as np
from scipy.integrate import solve_ivp

# Constants
avogadro_number = 6.022e23  # atoms/mol
density_steel = 7.9  # g/cm³ (approximate for stainless steel)
time_seconds = 360 * 24 * 3600  # 360 days in seconds

# Cladding composition (wt.%)
composition = {
    "Fe": 100 - (15.0 + 15.0 + 1.5 + 1.5 + 0.9 + 0.4 + 0.09),  # Balance is Iron
    "Cr": 15.0,
    "Ni": 15.0,
    "Mo": 1.5,
    "Mn": 1.5,
    "Si": 0.9,
    "Ti": 0.4,
    "C": 0.09,
}
b_ppm = 60  # Boron content in ppm
molar_masses = {
    "Fe": 55.487,
    "Cr": 51.996,
    "Ni": 58.690,
    "Mn": 54.938,
    "Si": 28.086,
    "Mo": 95.940,
    "C": 12.011,
    "Ti": 47.867,
}

# Cross-sections (cm²)
cross_sections = {
    "B10": 623e-24,  # 10B(n,α)
    "Ti": 0.23e-24,  # Ti(n,α)
    "Fe": 0.23e-24,  # Fe(n,α)
    "Cr": 0.20e-24,  # Cr(n,α)
    "Ni58": 4.2e-24,  # Ni(n,α)
}

# Isotopic abundance
isotopic_abundance = {
    "Ni58": 0.683,  # Abundance of 58Ni
}

# Node-specific neutron flux data
nodes = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
heights = np.array([42.5, 127.5, 212.5, 297.5, 382.5, 467.5, 552.5, 637.5, 722.5, 807.5])  # mm
peak_factors = np.array([0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.498])
peak_flux = 6.1e15  # n/cm²/s at the peak node

# Compute number densities (atoms/cm³)
total_mass = sum(composition[element] / molar_masses[element] for element in composition)
n_total = density_steel * avogadro_number / total_mass

number_densities = {element: (composition[element] / molar_masses[element]) * n_total for element in composition}
number_densities["B10"] = b_ppm * 1e-6 * number_densities["Fe"]  # Boron in ppm
number_densities["Ni58"] = isotopic_abundance["Ni58"] * number_densities["Ni"]

# Define the Bateman equations for helium production
def helium_production(t, y, flux, n_densities, cross_sections):
    n_he = y[0]  # Helium concentration
    production_rate = (
        cross_sections["B10"] * flux * n_densities["B10"]
        + cross_sections["Fe"] * flux * n_densities["Fe"]
        + cross_sections["Cr"] * flux * n_densities["Cr"]
        + cross_sections["Ti"] * flux * n_densities["Ti"]
        + cross_sections["Ni58"] * flux * n_densities["Ni58"]
    )
    return [production_rate]

# Solve helium concentration for each node
helium_results = []
for i, peak_factor in enumerate(peak_factors):
    flux = peak_flux * peak_factor  # Node-specific neutron flux
    y0 = [0]  # Initial helium concentration
    t_span = (0, time_seconds)  # Irradiation duration
    solution = solve_ivp(
        helium_production,
        t_span,
        y0,
        args=(flux, number_densities, cross_sections),
        method="RK45",
    )
    helium_results.append(solution.y[0][-1])  # Final helium concentration

# Convert helium concentrations to atomic ppm
helium_ppm = [he / n_total * 1e6 for he in helium_results]

# Calculate the weighted average of helium concentration
# Weighted by peak flux (or node height)
weights = peak_factors  # Using peak factor as weights for flux distribution
weighted_helium_avg = np.average(helium_ppm, weights=weights)

# Print results for each node
for i, ppm in enumerate(helium_ppm):
    print(f"Node {nodes[i]} (Height {heights[i]} mm): Helium concentration = {ppm:.2f} ppm")

# Print the weighted average helium concentration
print(f"Weighted average helium concentration: {weighted_helium_avg:.2f} ppm")
