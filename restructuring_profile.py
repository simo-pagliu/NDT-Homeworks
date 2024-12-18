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
plt.xlim(r_0, r_4)
plt.ylabel('Temperature [K]')

# Add legend to the plot
plt.legend()

# Put the legend out of the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show the figure
plt.show()




def youngs_modulus_function(temperature):
  y_m = 202.7 - 0.08167*temperature
  return y_m

def poisson_ratio_function(temperature):
  p_r = 0.277 + 6*10**-5*temperature
  return p_r

r = sp.symbols('r')
sigma_r = sp.symbols('sigma_r', cls=sp.Function)
alfa = sp.symbols('alfa', real = True, nonnegative = True)
E = sp.symbols('E', real = True, nonnegative = True)
ni = sp.symbols('ni')
T = sp.symbols('T', cls=sp.Function)
Ro =  3.31*10**-3
Ri = 3.22*10**-3
r_vector = np.linspace(Ri,Ro,50)
To = f.get_temperature_at_point(0.425, Ro, vars.T_map) - 273
Ti = f.get_temperature_at_point(0.425, Ri, vars.T_map) - 273
print(f'the value of temperature at inner cladding is {Ti}')
E_value1 = youngs_modulus_function(To)
E_value2 = youngs_modulus_function(Ti)
ni_value1 = poisson_ratio_function(To)
ni_value2 = poisson_ratio_function(Ti)

E_value_real = (E_value1 + E_value2)/2
ni_value_real = (ni_value1 + ni_value2)/2

#function definition, using superposition problem
func1 = sp.Eq(
    (r**3*sigma_r(r).diff(r)).diff(r),0
)

func2 = sp.Eq(
    ((alfa*E*r**2)/(1-ni))*T(r).diff(r),0
)
display(func1)
display(func2)

#C1,C2,C3,C4 = sp.symbols("C1 C2 C3 C4")
solution1 = sp.dsolve(func1,ics={sigma_r(r).subs(r, Ro): -0.2, sigma_r(r).subs(r, Ri): -2.5}) #pressure at the two extreme of the cladding
solution2 = sp.dsolve(func2,ics={T(r).subs(r, Ro): To}) #temperature condition at two cladding extreme
#ics={T(r).subs(r, Ro): To, T(r).subs(r, Ri): Ti}

#Print of the solutions
print("Solution of sigma_r(r):")
display(solution1)

print("\nSolution of T(r):")
display(solution2)

#Symbolic solution for sigma_theta
sigma_r_solution = solution1.rhs        #+ solution2.rhs
print("real solution of sigma_r:")
display(sigma_r_solution)
sigma_theta = sp.symbols('sigma_theta', cls=sp.Function)
func3 = sp.Eq(sigma_theta(r), sigma_r(r) + r*sigma_r(r).diff(r))
#func3 = sp.Eq(sigma_theta(r), sigma_r_solution + r * sp.diff(sigma_r_solution, r)) #l'errore sta qui nel calcolo della derivata di sigma_r
sigma_theta_eq = func3.rhs

print("Solution of sigma_theta:")
display(func3)

#Symbolic solution for sigma_z
sigma_z = sp.symbols('sigma_z', cls=sp.Function)
func4 = sp.Eq(sigma_z(r), sigma_r_solution + sigma_theta_eq)
sigma_z_eq = func4.rhs

print("Solution of sigma_z:")
display(func4)


#numerical solution for the plot
numerical_sigma_r = sp.lambdify(r, sigma_r_solution, 'numpy')
#numerical_sigma_theta = sp.lambdify(r, sigma_theta_eq, 'numpy')
#numerical_sigma_z = sp.lambdify(r,sigma_z_eq, 'numpy')

##Starting from sigma_r's equation, we are able to find in numerical way sigma_theta and sigma_z  
numerical_sigma_theta2 = lambda r: 40.3773694061598 - (0.000444569716950827)/(r**2) + (2*0.000444569716950827)/(r**2)
numerical_sigma_z2 = lambda r: 40.3773694061596 - (0.000444569716950827)/(r**2) + 40.3773694061596 - (0.000444569716950827)/(r**2) + (2*0.000444569716950827)/(r**2)

#Function to plot
sigma_r_plot = numerical_sigma_r(r_vector)
#sigma_theta_plot = numerical_sigma_theta(r_vector)
#sigma_z_plot = numerical_sigma_z(r_vector)
sigma_theta_plot2 = numerical_sigma_theta2(r_vector)
sigma_z_plot2 = numerical_sigma_z2(r_vector)

#first correction for the sigma_z's plot
if np.isscalar(sigma_z_plot):
    sigma_z_plot = np.full_like(r_vector, sigma_z_plot)

#plot of the graph
plt.plot(r_vector,sigma_r_plot, color='red', label = 'sigma_r'), plt.plot(r_vector, sigma_theta_plot2, color='green', label='sigma_theta'), plt.plot(r_vector, sigma_z_plot2, color='yellow',label = 'sigma_z')
plt.title("pipe equation's solution")
plt.xlabel('r (m)')
plt.ylabel('stress')
plt.grid(True)
plt.legend(loc='best')
plt.show()