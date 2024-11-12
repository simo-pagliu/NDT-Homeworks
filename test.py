# 3rd Party Libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import math
import dill
from IPython.display import display, Math
# Internal Libraries
import funs as f
import nuclei_func as nf

# %% LOAD DATA
# Load all objects in a single line using dill
with open("data_file.dill", "rb") as file:
    data_objects = dill.load(file)
# Access the objects with their original names
Cladding_Proprieties = data_objects["Cladding_Proprieties"]
Fuel_Proprieties = data_objects["Fuel_Proprieties"]
Coolant_Proprieties = data_objects["Coolant_Proprieties"]
Helium_Proprieties = data_objects["Helium_Proprieties"]
ThermoHydraulics = data_objects["ThermoHydraulics"]
Geometrical_Data_Cold = data_objects["Geometrical_Data"]

# %% Compute the burnup, just once since it's constant
# Burnup
# s --> days
# W/m * m --> W
# m^3
# g/cm^3 --> kg/m^3
MOX_volume = np.pi * (Geometrical_Data_Cold.fuel_outer_diameter[0]/2)**2 * Geometrical_Data_Cold.h_values[-1]
MOX_mass = MOX_volume* Fuel_Proprieties.Density*1000
MOX_molar_mass = nf.mixture(Fuel_Proprieties.Molar_Mass, Fuel_Proprieties.Qualities) * 1e-3 # g/mol --> kg/mol
MOX_mols = MOX_mass / MOX_molar_mass
HM_molar_mass= nf.mixture(Fuel_Proprieties.Molar_Mass.remove(2), Fuel_Proprieties.Qualities.remove(2))
HM_mass = MOX_mols * HM_molar_mass #kg of Heavy Metal

Burnup = ThermoHydraulics.uptime / (24*3600) \
        * ThermoHydraulics.q_linear_avg * Geometrical_Data_Cold.h_values[-1] \
        / HM_mass
Burnup = Burnup * 1e-6  # Wd/kgU * 1e-9/1e-3 --> GWd/t(HM) 
print("Bunrup", Burnup)

# %% MAIN LOOP
cladding_thicknesses = np.linspace(565e-6, 1e-6, 6)
T_fuel_out = 1160 # K (Initial guess)
iterations = []

# Loop over cladding thicknesses
for delta in cladding_thicknesses:
    # Initialize
    residual = 1  # Placeholder for residual
    j = 0
    previous_T_map = None  # Placeholder for previous T_map
    print("\n")
    # Loop until the residual is less than 1e-3
    while residual > 1e-3:
        j += 1
        Geometrical_Data = dill.loads(dill.dumps(data_objects["Geometrical_Data"]))  # Reset Geometrical_Data
        Geometrical_Data.thickness_cladding = np.full_like(Geometrical_Data.h_values, delta)
        print(f"Cladding Thickness: {Geometrical_Data.thickness_cladding[0]*1e3:.2f} mm - Iteration {j}", end='\r')
        # Verify if the fuel diameter is greater than the cladding diameter
        r_cladding_gap = np.array(Geometrical_Data.cladding_outer_diameter) / 2 - np.array(Geometrical_Data.thickness_cladding)
        r_gap_fuel = np.array(Geometrical_Data.fuel_outer_diameter) / 2
        diff = r_cladding_gap - r_gap_fuel

        if np.all(diff < 0):
            print(f"\033[91mFuel Diameter < Cladding Diameter, Gap is closed at iteration {j}\033[0m")
            break
        else:
            # Calculate T_map
            T_map = f.temperature_map(Coolant_Proprieties, Cladding_Proprieties, Helium_Proprieties, 
                                      Fuel_Proprieties, ThermoHydraulics, Geometrical_Data, T_fuel_out, Burnup)
            # We can do better with a vector of temperatures
            idx_fuel = np.argmin(np.abs(T_map.r[5, :] - Geometrical_Data.fuel_outer_diameter[0]/2))
            T_fuel_out = T_map.T[5, idx_fuel]
            # Thermal Expansion
            Geometrical_Data.fuel_outer_diameter, \
            Geometrical_Data.fuel_inner_diameter, \
            Geometrical_Data.cladding_outer_diameter,\
            Geometrical_Data.thickness_cladding \
            = f.thermal_expansion(Fuel_Proprieties, Cladding_Proprieties, Geometrical_Data_Cold, T_map)

            # Calculate residual if previous_T_map exists
            if previous_T_map is not None:
                residual = np.mean(np.abs(T_map.T - previous_T_map)) / np.mean(previous_T_map)
            previous_T_map = T_map.T.copy()
    iterations.append(T_map)

# Plot all iterations
plt.figure()
for T_map, delta in zip(iterations, cladding_thicknesses):
    plt.plot(T_map.r[5, :], T_map.T[5, :], label=f"{delta*1e3:.2f} mm")
plt.title("Temperature Profile")
plt.xlabel("Radius [m]")
plt.ylabel("Temperature [K]")
plt.legend()
plt.show()