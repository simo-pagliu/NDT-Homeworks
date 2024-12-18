# 3rd Party Libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import math
from IPython.display import display, Math

# Internal Libraries
import functions as f
import nuclei_func as nf
from functions import Material_Proprieties, ThermoHydraulicSpecs, GeometryData, DimensioningData, Temperature_Map



# %% [code]
#function to define the value of [wt%] in a precise radius
def get_plutonium_at_radius(radius_value):
  for i, D in enumerate(D_values):
    if D_values[i] == 0.1:
      Pu_r_value = Pu_r(radius_value,D)
  return Pu_r_value

# %% [code]
#coef definition #using the concept of conservation of plutonium
heights_of_slice_centre = [42.5, 127.5, 212.5, 297.5, 382.5, 467.5, 552.5, 637.5, 722.5, 807.5] # mm
D_values = [1,0.1,0.01,0.001,0.0001,0.00001] #exponential fit
Pu_start = 0.29 #reffered to the initial concentration (power) of plutonium
R = 2.71 #mm                                #(fuel_outer_diameter/2)*10**3
print(R)
R_void = 0.00038319
display(Math(r'R_{fuel} =' + f'{R:.6f}' + r'\text{mm}'))
def alfa_function(R_gas,E,T):
  alfa = E*10**3/(R_gas*T)
  return alfa
R_gas = 8.314 #J/molK
E = 350 #from paper(i have the name) it is written that the value is in the range of 350 and 450 kj/mol
T_ref = f.get_temperature_at_point(heights_of_slice_centre[4],0,vars.T_map) #starting the evaluation with the worst case
#alfa_value = alfa_function(R_gas,E,T_ref)
alfa_value = 10
display(Math(r'alfa =' + f'{alfa_value:.6f}' + r'\text{}'))
display(Math(r'T_{ref} =' + f'{T_ref:.6f}' + r'\text{ K}'))
radius_vector = np.linspace(R_void,R,1000)
r_star = R*0.207  #the star radius, which is the minimum of Pu concentration, is taken from the handouts

#Pu_r = lambda r,D: Pu_start + Pu_start*(D*(np.exp(-2*alfa_value*((r-r_star)/R)) + 2*np.exp(-alfa_value*((r-r_star)/R))))
#Pu_r = lambda r, D: 1 - D * (np.exp(-2 * alfa_value * ((r - r_star) / R)) - 2 * np.exp(-alfa_value * ((r - r_star) / R))) #ratio between the Pu concentration radius-dependent and the initial value of Pu concentration
Pu_r = lambda r,D: 1 + D * (+np.exp(-2 * alfa_value * (r - r_star) / R) - 2 * np.exp(-(alfa_value) * (r - r_star) / R))
#Pu_r = lambda r,D: 1 + D*(np.exp(-((r-r_star)**2)/R**2) - np.exp(-alfa_value*(r/R)))



# Generate and plot graphs for each value of D
for i, D in enumerate(D_values):
    Pu_concentration = Pu_r(radius_vector, D)
    plt.figure(figsize=(8, 5))
    plt.plot(radius_vector, Pu_concentration, label=f'D = {D}')
    plt.axvspan(0,R_void, color='lightcyan', alpha=0.3, label='Void region')
    plt.axhline(y=1, color='red', linestyle='--',linewidth = '3', label='start Pu concentration')
    plt.xlim(0, R)  # 0 as the lower bound, None for no upper bound

    # plotting of the graph
    plt.xlabel('Radius')
    plt.ylabel('Pu Concentration')
    plt.title(f'Plutonium Concentration vs Radius for D = {D}')
    plt.legend(loc = 'best')
    plt.grid(True)
    plt.show()