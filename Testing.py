import numpy as np
import nuclei_func as nf
import math
from functions import Material_Proprieties, GeometryData
geom_data = GeometryData(
    fuel_outer_diameter=5.42 * 1e-3,  # m - GIVEN
    fuel_inner_diameter=0.00 * 1e-3,  # m
    cladding_outer_diameter=6.55 * 1e-3,  # m - GIVEN
    thickness_cladding=0.26 * 1e-3, # m
    pin_pitch=8.275 * 1e-3,  # m
    h_values = np.linspace(0, 0.85, 1000), # m
    fuel_pellet_height = 7e-3, # m
    fuel_roughness = 0.5e-6, # m ---> TO BE COMPUTED PROPERLY
    cladding_roughness = 0.5e-6 # m ---> TO BE COMPUTED PROPERLY
)

mol_qual = nf.w2mol([0.711, 0.29], [235 + 2*16, 239 + 2*16])  # UO2, PuO2

fuel = Material_Proprieties(
    Elements=["U-235", "U-238", "O-16", "Pu"],
    Qualities=[mol_qual[0] * (1 - mol_qual[1]), (1 - mol_qual[0]) * (1 - mol_qual[1]), 2, mol_qual[1]], # Molar fractions
    Theoretical_Density=11.31, # g/cm^3
    Percent_of_Theoretical_Density = 94.5, # %
    Molar_Mass=[235.0439299, 238.05078826, 15.99491461956, 244.064204],  # g/mol
    Thermal_Conductivity=lambda k_inf, beta: 1.755 + (k_inf - 1.755) * math.exp(-beta),  # W/m K
    Emissivity = 0.79, # -
    Thermal_Expansion_Coeff=1.2e-5,  # 1/°C
    Specific_Heat=270,  # Approximate value in J/kg K for MOX fuel
    Melting_Temperature=lambda pu, x, beta: 2964.92 + ((3147 - 364.85 * pu - 1014.15 * x) - 2964.92) * math.exp(-beta),  # K
    Oxigen_to_metal_ratio = 1.957, # -
    Grain_diameter = 10 * 1e-6,  # m
    Youngs_Modulus=lambda t, p: (22.43 * 10**4 - 31.19 * t) * (1 - 2.6 * p),  # MPa
    Poissons_Ratio=0.32,  # Dimensionless
)

cladding = Material_Proprieties(
    Elements=["Cr", "Ni", "Mo", "Mn", "Si", "Ti", "C", "B [ppm]"],
    Qualities=[15, 15, 1.5, 1.5, 0.9, 0.4, 0.09, 60],
    Density=lambda eth: 7900 * (1 + eth)**-3,  # kg/m^3
    Thermal_Conductivity=lambda t: 13.95 + 0.01163 * t,  # W/m K
    Emissivity = 0.32, # -
    Thermal_Expansion_Coeff=lambda t: -3.101e-4 + 1.545e-5 * t + 2.75e-9 * t**2,  # 1/°C
    Specific_Heat=500,  # Approximate value in J/kg K for steel
    Melting_Temperature=1673,  # K
    
    Youngs_Modulus=lambda t: 202.7 - 0.08167 * t,  # GPa
    Poissons_Ratio=lambda t: 0.277 + 6e-5 * t,  # Dimensionless
    Yield_Stress=lambda t: 555.5 - 0.25 * t if t < 600 else (405.5 - 0.775 * (t - 600) if t < 1000 else 345.5 - 0.25 * t),  # MPa
    Ultimate_Tensile_Strength=lambda t: 700 - 0.3125 * t if t < 600 else (512.5 - 0.969 * (t - 600) if t < 1000 else 437.5 - 0.3125 * t),  # MPa
)

# Material: Helium (Filling Gas)
helium = Material_Proprieties(
    Elements=["He"],
    Qualities=[1],
    Density=0.1786,  # kg/m^3 at STP
    Thermal_Conductivity=lambda t: 15.8e-4 * t**0.79,  # W/m K
    Specific_Heat=5193,  # J/kg K at constant pressure
    Thermal_Expansion_Coeff=3.66e-3  # Approximate value for helium in 1/°C
)
############################################
# Functions
############################################
Burnup = 0
def k(temp):
    A = 0.01926 + 1.06e-6 * fuel.Oxigen_to_metal_ratio + 2.63e-8 * fuel.Molar_Mass[-1]
    B = 2.39e-4 + 1.37e-13 * fuel.Molar_Mass[-1]
    D = 5.27e9
    E = 17109.5

    # Calculate k_0
    k_0 = (1 / (A + B * temp) + D / (temp**2) * np.exp(-E / temp)) * (1 - fuel.Porosity)**2.5

    # Calculate final thermal conductivity k
    k = 1.755 + (k_0 - 1.755) * np.exp(-Burnup / 128.75)
    return k

# def thermal_resistance_gap(geom_data, helium, fuel, cladding, temperature_gas, temperature_fuel_out):
def k_gap(temperature_gas):
    # radius_gap_in = geom_data.fuel_outer_diameter / 2

    # Conductive heat transfer
    conduction = helium.Thermal_Conductivity(temperature_gas) / 0.5 # gap thickness
    temperature_fuel_out = 675 # !!!!!!!!!!!!!!!!!
    # Radiative heat transfer
    radiation = 4 * 5.67e-8 * (temperature_fuel_out**3) / (1/(fuel.Emissivity) + 1/(cladding.Emissivity) - 1)

    return conduction + radiation * 1 # gap outer radius

############################################
# Placeholder values
############################################
r_out_gap = 1
r_out_fuel = 0.5
T_r_out = 600
temp = 900
power = 20e3 # 20 kW

############################################
# Calculate temperature profile
############################################
temp_vals_gap = []
r_vect_gap = np.linspace(r_out_gap, r_out_fuel, 100)
for r in r_vect_gap:
    temp = T_r_out + power / (4 * np.pi * k_gap(temp)) * (r_out_gap**2 - r**2)
    temp_vals_gap.append(temp)

T_r_out = temp_vals_gap[99]

temp_vals_fuel = []
r_vect_fuel = np.linspace(r_out_fuel, 0, 100)
for r in np.linspace(r_out_fuel, 0, 100):
    temp = T_r_out + power / (4 * np.pi * k(temp)) * (r_out_fuel**2 - r**2)
    temp_vals_fuel.append(temp)


############################################
# Plotting
############################################
import matplotlib.pyplot as plt
plt.plot(r_vect_fuel, temp_vals_fuel, label="Fuel")
plt.plot(r_vect_gap, temp_vals_gap, label="Gap")
plt.legend()
plt.xlabel("Radius [m]")
plt.ylabel("Temperature [°C]")

plt.show()



####
    # for j, r in enumerate(r_plot[1:], start=1):       
    #     # In the pellet
    #     if r < geom_data.fuel_outer_diameter/2:
    #         th_res = thermal_resistance_fuel(Burnup, T_radial[j-1], fuel_data)
    #         r_out = geom_data.fuel_outer_diameter / 2
    #         # Compute the temperature
    #         T_value = T_radial[idx_fuel_r] + power * th_res * (r_out**2 - r**2)
        
    #     # In the gap
    #     elif r < geom_data.cladding_outer_diameter / 2 - geom_data.thickness_cladding:
    #         r_out = geom_data.cladding_outer_diameter / 2 - geom_data.thickness_cladding
    #         th_res = thermal_resistance_gap(geom_data, helium_data, fuel_data, cladding_data, T_radial[j-1], T_radial[idx_gap_r-1])         
    #         # Compute the temperature
    #         T_value = T_radial[idx_gap_r] + power * th_res * (r_out**2 - r**2)

    #     # Compute new value of T
    #     T_radial.append(T_value)