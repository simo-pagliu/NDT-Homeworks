import numpy as np
from math import log
import matplotlib.pyplot as plt
import sympy as sp

import functions as f
import nuclei_func as nf
from functions import Material_Proprieties, ThermoHydraulicSpecs, GeometryData, DimensioningData, Temperature_Map
from IPython.display import display, Math

#  2D plot of the temperature profile


R = GeometryData.fuel_outer_diameter/2
R_void = 0.511873*10**-3
h_vals = GeometryData.h_values
R_new, R_start, T_hot = f.cold_to_hot_fuel(Material_Proprieties,GeometryData,vars,h_vals) #fuel R after thermal expansion
R_new_value = R_new[0]
print(R_new_value)
r_fuel_vector = np.linspace(R_void,R_new,1000)

r_gap_fuel = R_new_value
r_end = R_void

#x_axis definition
r_gap_fuel = GeometryData.fuel_outer_diameter / 2
r_end = GeometryData.fuel_inner_diameter / 2
r_plot = np.linspace(r_gap_fuel, r_end, 25)
r_fuel = np.linspace(R_void,R_new_value)

#plot fter restructuring
temp_plot_bottom2 = [(f.get_temperature_at_point(0, r, vars.T_map)*(1-((2*R_void**2)/(R_new_value**2-R_void**2))*log(R_new_value/R_void))) for r in r_fuel]
temp_plot_center2 = [(f.get_temperature_at_point(0.425, r, vars.T_map)*(1-((2*R_void**2)/(R_new_value**2-R_void**2))*log(R_new_value/R_void))) for r in r_fuel]
temp_plot_top2 = [(f.get_temperature_at_point(0.850, r, vars.T_map)*(1-((2*R_void**2)/(R_new_value**2-R_void**2))*log(R_new_value/R_void))) for r in r_fuel]

# 2D plot of the temperature profile before restructuring
temp_plot_bottom1 = [f.get_temperature_at_point(0, r, vars.T_map) for r in r_plot]
temp_plot_center1 = [f.get_temperature_at_point(0.425, r, vars.T_map) for r in r_plot]
temp_plot_top1 = [f.get_temperature_at_point(0.850, r, vars.T_map) for r in r_plot]

# Create plot
plt.plot(r_plot * 1e3, temp_plot_bottom2, label='Bottom w\ restructuring', marker='o')
plt.plot(r_plot * 1e3, temp_plot_center2, label='Center w\ restructuring', marker='o')
plt.plot(r_plot * 1e3, temp_plot_top2, label='Top w\ restructuring', marker='o')

# Create plot
plt.plot(r_plot*10**3 , temp_plot_bottom1, label='Bottom w\o restructuring', marker='o', color ='black')
plt.plot(r_plot*10**3 , temp_plot_center1, label='Center w\o restructuring', marker='o', color ='grey')
plt.plot(r_plot*10**3 , temp_plot_top1, label='Top w\o restructuring', marker='o', color ='brown')

r_0 = GeometryData.fuel_inner_diameter/2 * 1e3
r_1 = GeometryData.fuel_outer_diameter/2 * 1e3

# Add shading to different regions
colors = ['#00008B', '#0000CD', '#4169E1', '#6495ED', '#87CEEB']
plt.axvspan(0,R_void*10**3, color='blue', alpha=0.3, label='Void region')
plt.axvspan(r_0, r_1, color=colors[4], alpha=0.3, label='Fuel Region')


# Set title and axis labels
plt.title('Temperature Profile')
plt.xlabel('Radius [mm]')
plt.xlim(r_0, r_1)
plt.ylabel('Temperature [K]')

# Add legend to the plot
plt.legend()

# Put the legend out of the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show the figure
plt.show()