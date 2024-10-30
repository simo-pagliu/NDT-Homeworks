# 3rd Party Libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import math
import dill
from IPython.display import display, Math

# Internal Libraries
import functions_looping as f
import nuclei_func as nf
from functions_looping import Material_Proprieties, ThermoHydraulicSpecs, GeometryData, DimensioningData, Temperature_Map

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
fig.show()
