import sympy as sp
import numpy as np
from scipy.optimize import fsolve

def void_geometry(R_fuel,R_equiaxed,R_columnar,density_equiaxed_ratio,density_columnar_ratio, density_TD):
    #the main problem is the determination of equiaxed and columnar radius

    R_void = sp.symbols('R_void')
    R_void = sp.solvers.solve(R_fuel**2*density_TD-(R_columnar**2-R_void**2)*(density_columnar_ratio*density_TD)-(R_equiaxed**2-R_columnar**2)*(
            density_equiaxed_ratio*density_TD)-(R_fuel**2-R_equiaxed**2)*(density_TD))[0]

    return R_void 

# Function to calculate A based on x and Pu_concentration
def A_f(x, Pu_concentration):
    A_value = 0.01926 + 1.06 * 10**-6 * x + 2.63 * 10**-8 * Pu_concentration
    return A_value

# Function to calculate B based on Pu_concentration
def B_f(Pu_concentration):
    B_value = 2.39 * 10**-4 + 1.37 * 10**-13 * Pu_concentration
    return B_value

# Function to calculate porosity based on density
def porosity_f(density_value):
    p = 1 - density_value  # porosity evaluation in different zones
    return p

# Burnup calculation
def bu_function(fuel_outer_diameter, fuel_height, fuel_density, uptime, q_linear_avg):
  Bu = uptime / (24 * 3600) * q_linear_avg * fuel_height / (np.pi * (fuel_outer_diameter / 2)**2 * fuel_height * fuel_density)
  Bu = Bu * 1e-6  # Conversion to GWd/tHM
  return Bu

# Definition of symbols
T = sp.symbols('T')
A = sp.symbols('A')
B = sp.symbols('B')
D = sp.symbols('D')
E = sp.symbols('E', positive=True)
p = sp.symbols('p')
k = sp.Function('k')(T)

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
D_value = 5.27 * 10**9
E_value = 17109.5
# Data needed for burnup calculation
fuel_outer_diameter = 5.42e-3  # m
fuel_height = 0.85  # m (total fuel height)
fuel_density = 11.31 * 0.945 * 1000  # kg/m^3 (theoretical density in g/cm^3 * percentage of theoretical density * conversion to kg/m^3)
uptime = 360 * 24 * 3600  # s (operating time in seconds)
q_linear_avg = 38.7e3  # W/m (average linear power)
# Temperatures
T3 = 1400  # maximum temperature at the outer fuel radius
T2 = 1600 + 273  # temperature at the equiaxed radius
T1 = 1800 + 273  # temperature at the columnar radius

Burnup = bu_function(fuel_outer_diameter, fuel_height, fuel_density, uptime, q_linear_avg)

# Print result
print("Burnup:", Burnup, "GWd/tHM")

# Equation to integrate to obtain thermal conductivity
k_function = sp.Eq(k, 1.755 + (((1 / (A + B * T) + (D / T**2) * sp.exp(-E / T) * (1 - p)**2.5)) - 1.755) * sp.exp(-Burnup / 128.75)) # the problem in this equation is the burnup
Equation_to_integrate = k_function.rhs

# Symbolic integration
integral_k = sp.integrate(Equation_to_integrate, T)

# Creation of the matrix for the integrated values of k
n = 1
m = 3
matrix_of_integrals = sp.MutableDenseMatrix(n, m, [0] * n * m)

# Print the symbolic integral of k(T)
print("Symbolic integral of k(T) is:")
sp.pprint(integral_k)

# Calculation and substitution of values in the matrix cells
for i in range(m):
    matrix_of_integrals[0, i] = integral_k.subs({
        A: A_value,
        B: B_value,
        D: D_value,
        E: E_value,
        p: porosity_f(DENSITY[i])
    })

# Print the matrix
print(matrix_of_integrals)

matrix_of_integrals_numeric = sp.lambdify(T, matrix_of_integrals, 'numpy')
print(f'The matrix is: {matrix_of_integrals_numeric}')

# Convert the symbolic integral to a numerical function for evaluation with temperature values

k3_start = matrix_of_integrals[0, 0].subs({T: T3})
k3_end = matrix_of_integrals[0, 0].subs({T: T2})
k2_start = matrix_of_integrals[0, 1].subs({T: T2})
k2_end = matrix_of_integrals[0, 1].subs({T: T1})
integral_k_numeric = sp.lambdify(T, matrix_of_integrals, 'numpy')

# Calculate numerical values for k3_start, k3_end, k2_start, k2_end
k3_start = integral_k_numeric(T3)[0, 0]
k3_end = integral_k_numeric(T2)[0, 0]
k2_start = integral_k_numeric(T2)[0, 1]
k2_end = integral_k_numeric(T1)[0, 1]

print("k3_start:", k3_start)
print("k3_end:", k3_end)
print("k2_start:", k2_start)
print("k2_end:", k2_end)

# System of equations
R_fuel = 0.542/2
q_first_avg = q_linear_avg
print(f'the value of q linear is{q_first_avg}')

def system(variables):
    R_void, R_columnar, R_equiaxed = variables

    # Equation 1
    eq1 = (np.pi * R_fuel**2 * Fuel_density -
           (np.pi * (R_columnar**2 - R_void**2) * columnar_density +
            np.pi * (R_equiaxed**2 - R_columnar**2) * equiaxed_density +
            np.pi * (R_fuel**2 - R_equiaxed**2) * 0.88 * Fuel_density))

    # Equation 2
    eq2 = (k3_end - k3_start - (q_first_avg / (4 * np.pi)) * (1 - (R_equiaxed / R_fuel)**2))

    # Equation 3
    eq3 = (k2_end - k2_start - 
           (q_first_avg / (4 * np.pi)) * (columnar_density / (0.88 * density_TD)) *
           ((R_columnar**2 - R_void**2) / R_fuel**2) *
           (1 - (R_columnar / R_equiaxed)**2 - ((equiaxed_density - 0.88 * Fuel_density) / equiaxed_density) *
            np.log(R_equiaxed / R_columnar)**2))

    return [eq1, eq2, eq3]

# Initial values for variables
initial_values = [0.17, 0.35, 0.47]   #unica problema sta qui con i valori iniziali, che tipicamente sono quelli che poi stampa a terminale come risultato

# Solving the system of equations
numerical_solutions = fsolve(system, initial_values)

#this is a bond on the result
if all(0 <= sol <= R_fuel for sol in numerical_solutions) and \
   numerical_solutions[0] < numerical_solutions[1] < numerical_solutions[2] < R_fuel:
    print("The solutions are consistent with the image constraints.")
else:
    print("Warning: the solutions are not consistent with the maximum size constraints.")

# Print the solutions (Numerical way to find POSSIBLE SOLUTIONS for the problem)
print("Numerical solutions:")
print("R_void:", numerical_solutions[0])
print("R_columnar:", numerical_solutions[1])
print("R_equiaxed:", numerical_solutions[2])

#def void_geometry(R_fuel,R_equiaxed,R_columnar,density_equiaxed_ratio,density_columnar_ratio.density_TD)
#the main problem is the determination of equiaxed and columnar radius



















#mass equation
equation1 = sp.Eq(sp.pi*R_fuel**2*density_TD,sp.pi*(R_columnar**2-R_void**2)*columnar_density
                  + sp.pi*(R_equiaxed**2-R_columnar**2)*equiaxed_density +
                  sp.pi*(R_fuel**2-R_equiaxed**2)*0.88*density_TD
)

#heat conduction equations   #è ancora una bozza, inizio a inserirla nel codice
k3_start = matrix_of_integrals[0,0].subs({T:T3})
k3_end = matrix_of_integrals[0,0].subs({T:T2})
equation2 = sp.Eq(k3_end - k3_start,(q_first_avg/(4*sp.pi))*(1-(R_equiaxed/R_fuel)**2))

k2_start = matrix_of_integrals[0,1].subs({T:T2})
k2_end = matrix_of_integrals[0,1].subs({T:T1})
equation3 = sp.Eq(k3_end - k3_start,(q_first_avg/(4*sp.pi))*(columnar_density/(0.88*density_TD))*((R_columnar**2-R_void**2)/R_fuel**2)*
                  (1-(R_columnar/R_equiaxed)**2 - ((equiaxed_density-0.88*density_TD)/equiaxed_density)*sp.log(R_equiaxed/R_columnar)**2))
solution = sp.solvers.solve((equation1, equation2, equation3), (R_void, R_columnar, R_equiaxed))
print(solution)