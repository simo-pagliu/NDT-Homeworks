import dill
import math
import nuclei_func as nf
import numpy as np
from loop import Material_Proprieties, ThermoHydraulicSpecs, GeometryData
# Material: Cladding
# 15-15, Ti stabilized, cold worked stainless steel
Cladding_Proprieties = Material_Proprieties(
    Elements=["Cr", "Ni", "Mo", "Mn", "Si", "Ti", "C", "B [ppm]"],
    Qualities=[15, 15, 1.5, 1.5, 0.9, 0.4, 0.09, 60],
    Density=lambda eth: 7900 * (1 + eth)**-3,  # kg/m^3
    Thermal_Conductivity=lambda t: 13.95 + 0.01163 * t,  # W/m K
    Emissivity = 0.32, # -
    Thermal_Expansion_Coeff=lambda t: -3.101e-4 + 1.545e-5 * (t-273.15) + 2.75e-9 * (t-273.15)**2,  # 1/°C --> 1/K
    Specific_Heat=500,  # Approximate value in J/kg K for steel
    Melting_Temperature=1673,  # K
    
    Youngs_Modulus=lambda t: 202.7 - 0.08167 * t,  # GPa
    Poissons_Ratio=lambda t: 0.277 + 6e-5 * t,  # Dimensionless
    Yield_Stress=lambda t: 555.5 - 0.25 * t if t < 600 else (405.5 - 0.775 * (t - 600) if t < 1000 else 345.5 - 0.25 * t),  # MPa
    Ultimate_Tensile_Strength=lambda t: 700 - 0.3125 * t if t < 600 else (512.5 - 0.969 * (t - 600) if t < 1000 else 437.5 - 0.3125 * t),  # MPa
)

# Material: Fuel
# Homogeneous MOX fuel

# Qualities have to be converted from weight to molar (won't change in hot condition):
mol_qual = nf.w2mol([0.711, 0.29], [235 + 2*16, 239 + 2*16])  # UO2, PuO2

Fuel_Proprieties = Material_Proprieties(
    Fission_Yield=0.3,
    Elements=["U-235", "U-238", "O-16", "Pu"],
    Qualities=[mol_qual[0] * (1 - mol_qual[1]), (1 - mol_qual[0]) * (1 - mol_qual[1]), 2, mol_qual[1]], # Molar fractions
    Micro_Fission = [1.047756375, 0.55801001, 0, 1.689844625],  # barn
    Theoretical_Density=11.31, # g/cm^3
    Percent_of_Theoretical_Density = 94.5, # %
    Porosity_Columnar = 0.02,
    Porosity_Equiaxed = 0.04,
    Molar_Mass=[235, 238, 16, 239],  # g/mol
    Thermal_Conductivity=lambda k_inf, beta: 1.755 + (k_inf - 1.755) * math.exp(-beta),  # W/m K
    Emissivity = 0.79, # -
    Thermal_Expansion_Coeff=1.2e-5,  # 1/°C
    Specific_Heat=270,  # Approximate value in J/kg K for MOX fuel
    Melting_Temperature=lambda pu, x, beta: 2964.92 + ((3147 - 364.85 * pu - 1014.15 * x) - 2964.92) * math.exp(-beta/41.01),  # K
    Oxigen_to_metal_ratio = 1.957, # -
    Grain_diameter = 10 * 1e-6,  # m
    Youngs_Modulus=lambda t, p: (22.43 * 10**4 - 31.19 * t) * (1 - 2.6 * p),  # MPa
    Poissons_Ratio=0.32,  # Dimensionless
)

# Material: Coolant (Sodium)
Coolant_Proprieties = Material_Proprieties(
    Elements=["Na"],
    Qualities=[1],
    Density=lambda t: 954.1579 + ((t-273) * 9/5 +32) *( ((t-273) * 9/5 +32) * ((((t-273) * 9/5 +32) * 0.9667e-9 - 0.46e-5)) - 0.1273534),  # kg/m^3 (t * 9/5 +32) is the convertion from C to K to F
    Viscosity=lambda t: (np.exp(813.9 / t -2.530 ))/1000,  # Pa s
    Thermal_Conductivity=lambda t: 110 - 0.0648 * t + 1.16e-5 * t**2,  # W/m K
    Specific_Heat=lambda t: 1608 - 0.7481 * t + 3.929e-4 * t**2,  # J/kg K
    
    Melting_Temperature=98,  # °C
    Boiling_Temperature=882,  # °C
    Nusselt_Number=lambda pe: 7 + 0.025 * pe**0.8
)

# Material: Helium (Filling Gas)
Helium_Proprieties = Material_Proprieties(
    Elements=["He"],
    Qualities=[1],
    Density=0.1786,  # kg/m^3 at STP
    Thermal_Conductivity=lambda t, x: (15.8e-4 * t**0.79)**x * (0.935e-4 * t**0.79)**(1-x),  # W/m K
    Specific_Heat=5193,  # J/kg K at constant pressure
    Thermal_Expansion_Coeff=3.66e-3  # Approximate value for helium in 1/°C
)

len_h = 10
Geometrical_Data = GeometryData(
    Initial_Gas_Temperature=20 + 273,  # K
    Initial_Gas_Pressure=1e5,  # Pa
    fuel_outer_diameter=[5.42 * 1e-3]*len_h,  # m - GIVEN
    fuel_inner_diameter=[0.00 * 1e-3]*len_h,  # m
    cladding_outer_diameter=[6.55 * 1e-3]*len_h,  # m - GIVEN
    thickness_cladding=[0.3 * 1e-3]*len_h, # m
    pin_pitch=8.275 * 1e-3,  # m
    h_values = np.linspace(0, 0.85, len_h), # m
    fuel_pellet_height = 7e-3, # m
    fuel_roughness = 2e-6, # m
    cladding_roughness = 1e-6 # m
)

# Example of initializing Thermo-Hydraulic specifications
heights_of_slice_centre = [42.5, 127.5, 212.5, 297.5, 382.5, 467.5, 552.5, 637.5, 722.5, 807.5] # mm
ThermoHydraulics = ThermoHydraulicSpecs(
    coolant_inlet_temp=395 + 273,  # K
    coolant_inlet_pressure=1e5,  # Pa
    coolant_mass_flow_rate=0.049,  # kg/s
    q_linear_avg = 38.7e3,  #W/m,
    uptime = 360 * 24 * 3600,  # s
    h_peak_factor = [h * 1e-3 for h in heights_of_slice_centre],  # m
    peak_factors = [0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.498],
    neutron_flux_peak = 6.1e15  # Neutron Flux (> 100 keV) (10^15 n cm^-2 s^-1) at Peak Power Node    
)

data_objects = {
    "Cladding_Proprieties": Cladding_Proprieties,
    "Fuel_Proprieties": Fuel_Proprieties,
    "Coolant_Proprieties": Coolant_Proprieties,
    "Helium_Proprieties": Helium_Proprieties,
    "Geometrical_Data": Geometrical_Data,
    "ThermoHydraulics": ThermoHydraulics,
}

# Save all objects in a single file using dill
with open("data_file.dill", "wb") as file:
    dill.dump(data_objects, file)
