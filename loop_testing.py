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
from funs import DimensioningData

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
Geometrical_Data = dill.loads(dill.dumps(data_objects["Geometrical_Data"]))
Geometrical_Data_Cold = dill.loads(dill.dumps(data_objects["Geometrical_Data"]))

# %% Da sistemareeeee
temps = [ThermoHydraulics.coolant_inlet_temp]
T_fuel_out = 1160 # K (Initial guess)
vars = DimensioningData(
    filling_gas_pressure = 1e5,  # Pa
    filling_gas_temperature = 20,  # °C
    temperature_map = ''
)
# %% Compute the burnup, just once since it's constant
# Burnup
# s --> days
# W/m * m --> W
# m^3
# g/cm^3 --> kg/m^3
Burnup = ThermoHydraulics.uptime / (24*3600) \
        * ThermoHydraulics.q_linear_avg * Geometrical_Data.h_values[-1] \
        / (np.pi * (Geometrical_Data.fuel_outer_diameter[0]/2)**2 * Geometrical_Data.h_values[-1] \
        * Fuel_Proprieties.Density*1000)
Burnup = Burnup * 1e-6  # Wd/kgU * 1e-9/1e-3 --> GWd/t(HM) 
print("Bunrup", Burnup)

# %% Tempearture Map
# import time
# tic = time.time()
# toc = time.time()
# print(f"Elapsed time: {toc - tic} s")

# %% Plot meshgrid to chek it
# X, Y = f.create_meshgrid(Geometrical_Data)
# plt.figure()
# plt.plot(X, Y, 'o', markersize=0.5)
# plt.title('Meshgrid')
# plt.xlabel('Radius [m]')
# plt.ylabel('Height [m]')
# plt.show()

# %% Adding Thermal Expansion
cladding_thicknesses = np.linspace(565e-6, 1e-6, 6)
cladding_thicknesses = [100e-6]
plt.figure()
T_fuel_out_vector = []

for delta in cladding_thicknesses:
    Geometrical_Data.thickness_cladding = [delta]*len(Geometrical_Data.h_values)
    print(f"Cladding Thickness: {Geometrical_Data.thickness_cladding[0]*1e3:.2f} mm")
    for j in range(10):
        # Verify if the fuel diameter is greater than the cladding diameter
        r_cladding_gap = np.array(Geometrical_Data.cladding_outer_diameter) / 2 - np.array(Geometrical_Data.thickness_cladding)
        r_gap_fuel = np.array(Geometrical_Data.fuel_outer_diameter) / 2
        diff = r_cladding_gap - r_gap_fuel
        if np.all(diff < 0):
            print(f"\033[91mFuel Diameter < Cladding Diameter, Gap is closed at iteration {j}\033[0m")
            break
        else:
            # Execute core functions
        ############################################################################################################
        # Thermal map computation
            T_map = f.temperature_map(Coolant_Proprieties, Cladding_Proprieties, Helium_Proprieties, Fuel_Proprieties, \
                                            ThermoHydraulics, Geometrical_Data, T_fuel_out, Burnup)
            # We can do better with a vector of temperatures
            idx_fuel = np.argmin(np.abs(T_map.r[5, :] - Geometrical_Data.fuel_outer_diameter[0]/2))
            T_fuel_out = T_map.T[5, idx_fuel]
            T_fuel_out_vector.append(T_fuel_out)
        # Thermal Expansion
            Geometrical_Data.fuel_outer_diameter, \
            Geometrical_Data.fuel_inner_diameter, \
            Geometrical_Data.cladding_outer_diameter,\
            Geometrical_Data.thickness_cladding \
            = f.thermal_expansion(Fuel_Proprieties, Cladding_Proprieties, Geometrical_Data_Cold, T_map)
        ############################################################################################################

            # Plot profile
            T_this = T_map.T[5, :]
            plt.plot(T_map.r[5, :], T_this, '--', label=f"Cladding Thickness: {delta*1e3:.2f} mm at iteration {j}")

# Set title and axis labels
# plt.plot(T_map.r[5, :], T_this, '--', label=f"Cladding Thickness: {delta*1e3:.2f} mm")
plt.title('Temperature Profile')
plt.xlabel('Radius [mm]')
plt.ylabel('Temperature [K]')
# Add legend to the plot
plt.legend()
# Show the figure
plt.show()

# %% Plot Fuel Temperature
plt.figure()
plt.plot(T_fuel_out_vector)
plt.title('Fuel Temperature')
plt.xlabel('Iterations')
plt.ylabel('Temperature [K]')
plt.show()

# %% 3D Plot Temperature Map
# fig = go.Figure(data=[go.Surface(z=T_map.T, x=T_map.r*1e3, y=T_map.h*1e3, colorscale='Thermal')])

# # Update the layout
# fig.update_layout(
#     scene=dict(
#         xaxis_title='Radius (mm)',
#         yaxis_title='Height (mm)',
#         zaxis_title='Temperature (K)'
#     ),
#     title="3D Temperature Profile vs. Radius and Height (q_values)",
#     autosize=False,
#     width=800,
#     height=800
# )

# Show the plot
# fig.show()
