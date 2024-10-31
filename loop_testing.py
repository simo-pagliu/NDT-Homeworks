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
Geometrical_Data = data_objects["Geometrical_Data"]
ThermoHydraulics = data_objects["ThermoHydraulics"]

# %% Da sistemareeeee
temps = [ThermoHydraulics.coolant_inlet_temp]
Burnup = 360 * 38.7 * 0.85 / (np.pi * (Geometrical_Data.fuel_outer_diameter/2)**2 * 0.85 * 0.945 * 11310) * 1e-3  # kWd/kgU
T_fuel_out = 1500 # K (Initial guess)
vars = DimensioningData(
    filling_gas_pressure = 1e5,  # Pa
    filling_gas_temperature = 20,  # °C
    temperature_map = ''
)

# %% Tempearture Map
import time
tic = time.time()
T_map = f.temperature_map(Coolant_Proprieties, Cladding_Proprieties, Helium_Proprieties, Fuel_Proprieties, \
                                    ThermoHydraulics, Geometrical_Data, T_fuel_out, Burnup)
toc = time.time()
print(f"Elapsed time: {toc - tic} s")

# %% Loop over cladding thicknesses
cladding_thicknesses = np.linspace(565e-6, 1e-6, 6)
plt.figure()

for delta in cladding_thicknesses:
    Geometrical_Data.thickness_cladding = delta
    print(f"Cladding Thickness: {Geometrical_Data.thickness_cladding*1e3:.2f} mm")
    T_map = f.temperature_map(Coolant_Proprieties, Cladding_Proprieties, Helium_Proprieties, Fuel_Proprieties, \
                                    ThermoHydraulics, Geometrical_Data, T_fuel_out, Burnup)
    T_this = T_map.T[5, :]
    plt.plot(T_map.r[5, :], T_this, '--', label=f"Cladding Thickness: {delta*1e3:.2f} mm")

# Set title and axis labels
plt.title('Temperature Profile')
plt.xlabel('Radius [mm]')
plt.ylabel('Temperature [K]')
# Add legend to the plot
plt.legend()
# Show the figure
plt.show()

# %% 3D Plot Temperature Map
fig = go.Figure(data=[go.Surface(z=T_map.T, x=T_map.r*1e3, y=T_map.h*1e3, colorscale='Thermal')])

# Update the layout
fig.update_layout(
    scene=dict(
        xaxis_title='Radius (mm)',
        yaxis_title='Height (mm)',
        zaxis_title='Temperature (K)'
    ),
    title="3D Temperature Profile vs. Radius and Height (q_values)",
    autosize=False,
    width=800,
    height=800
)

# Show the plot
# fig.show()
