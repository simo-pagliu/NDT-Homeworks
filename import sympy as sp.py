import sympy as sp
import numpy as np
from scipy.optimize import fsolve
from IPython.display import display, Math

# Updated void_geometry function
def void_geometry(R_fuel, R_equiaxed, R_columnar, density_equiaxed_ratio, density_columnar_ratio, density_TD):
    # Determine R_void using symbolic solving
    R_void = sp.symbols('R_void', positive=True)
    eq = sp.Eq(R_fuel**2 * density_TD, 
               (R_columnar**2 - R_void**2) * (density_columnar_ratio * density_TD) + 
               (R_equiaxed**2 - R_columnar**2) * (density_equiaxed_ratio * density_TD) +
               (R_fuel**2 - R_equiaxed**2) * density_TD)
    R_void_solution = sp.solve(eq, R_void)
    
    # Ensure we return a valid solution (positive value)
    R_void_value = [sol.evalf() for sol in R_void_solution if sol.is_real and sol > 0]
    
    return R_void_value[0] if R_void_value else None

# Function to calculate A based on x and Pu_concentration
def A_f(x, Pu_concentration):
    A_value = 0.01926 + 1.06e-6 * x + 2.63e-8 * Pu_concentration
    return A_value

# Function to calculate B based on Pu_concentration
def B_f(Pu_concentration):
    B_value = 2.39e-4 + 1.37e-13 * Pu_concentration
    return B_value

# Function to calculate porosity based on density
def porosity_f(density_value, density_TD):
    p = 1 - (density_value / density_TD)  # Corrected porosity evaluation
    return p

# Burnup calculation
def bu_function(fuel_outer_diameter, fuel_height, fuel_density, uptime, q_linear_avg):
    Bu = uptime / (24 * 3600) * q_linear_avg * fuel_height / (np.pi * (fuel_outer_diameter / 2)**2 * fuel_height * fuel_density)
    Bu = Bu * 1e-6  # Conversion to GWd/tHM
    return Bu

# Initial data
density_TD = 11.31  # g/cm^3
Fuel_density = 11.31 * 0.945
equiaxed_density = 0.95 * density_TD
columnar_density = 0.98 * density_TD
DENSITY = [columnar_density, equiaxed_density, Fuel_density]
Pu_concentration = 0.29
O_M = 1.957
x = 2 - O_M
A_value = A_f(x, Pu_concentration)
B_value = B_f(Pu_concentration)
D_value = 5.27e9
E_value = 17109.5

# Data needed for burnup calculation
fuel_outer_diameter = 5.42e-3  # m
fuel_height = 0.85  # m
fuel_density = 11.31 * 0.945 * 1000  # kg/m^3
uptime = 360 * 24 * 3600  # s
q_linear_avg = 38.7e3  # W/m

# Burnup calculation
Burnup = bu_function(fuel_outer_diameter, fuel_height, fuel_density, uptime, q_linear_avg)

# Definition of symbols for thermal conductivity calculation
T = sp.symbols('T')
A = sp.symbols('A')
B = sp.symbols('B')
D = sp.symbols('D')
E = sp.symbols('E', positive=True)
p = sp.symbols('p')
k = sp.Function('k')(T)

# Equation to integrate to obtain thermal conductivity
k_function = sp.Eq(k, 1.755 + (((1 / (A + B * T) + (D / T**2) * sp.exp(-E / T) * (1 - p)**2.5)) - 1.755) * sp.exp(-Burnup / 128.75))
Equation_to_integrate = k_function.rhs

# Symbolic integration
integral_k = sp.integrate(Equation_to_integrate, T)

# Calculation and substitution of values in the matrix cells
moitaiav = []
for rho in DENSITY:
    moitaiav.append(integral_k.subs({
        A: A_value,
        B: B_value,
        D: D_value,
        E: E_value,
        p: porosity_f(rho, density_TD)
    }))

# Displaying results
display(Math(sp.latex(moitaiav[0])))
display(Math(sp.latex(moitaiav[1])))
display(Math(sp.latex(moitaiav[2])))

# Integral calculations
T_f_out = 1600  # max temperature at outer fuel radius
T_equiaxed = 1600 + 273  # temperature at equiaxed radius
T_columnar = 1800 + 273  # temperature at columnar radius

integral_as_fabricated = moitaiav[2].subs(T, T_f_out) - moitaiav[2].subs(T, T_equiaxed)
integral_equiaxed = moitaiav[1].subs(T, T_columnar) - moitaiav[1].subs(T, T_equiaxed)

display(Math(r"\text{Integral as fabricated }" + f"{integral_as_fabricated}"))
display(Math(r"\text{Integral equiaxed }" + f"{integral_equiaxed}"))

# Solve for R_equiaxed using equation 2
R_fuel = 0.542 / 2
q_first_avg = 38.7e3

eq2 = sp.Eq(integral_as_fabricated, (q_first_avg / (4 * np.pi)) * (1 - (sp.Symbol('R_equiaxed', positive=True) / R_fuel)**2))
R_equiaxed_solution = sp.solve(eq2, sp.Symbol('R_equiaxed', positive=True))
R_equiaxed_value = [sol.evalf() for sol in R_equiaxed_solution if sol.is_real and sol > 0][0]
display(Math(r"R_{equiaxed} = " + f"{R_equiaxed_value}"))

# Solve for R_void and R_columnar
R_columnar = sp.Symbol('R_columnar', positive=True)
R_void_value = void_geometry(R_fuel, R_equiaxed_value, R_columnar, 0.95, 0.98, density_TD)

# Solve for R_columnar using equation 3
eq3 = sp.Eq(integral_equiaxed,
            (q_first_avg / (4 * np.pi)) * (columnar_density / density_TD) *
            ((R_columnar**2 - R_void_value**2) / R_fuel**2) *
            (1 - (R_columnar / R_equiaxed_value)**2 - ((equiaxed_density - Fuel_density) / equiaxed_density) *
             (R_equiaxed_value / R_columnar - 1)**2))

R_columnar_solution = sp.solve(eq3, R_columnar)
R_columnar_value = [sol.evalf() for sol in R_columnar_solution if sol.is_real and sol > 0][0]
display(Math(r"R_{columnar} = " + f"{R_columnar_value}"))

R_void_value_final = R_void_value.subs(R_columnar, R_columnar_value)
display(Math(r"R_{void} = " + f"{R_void_value_final}"))
