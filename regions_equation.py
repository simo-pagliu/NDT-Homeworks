import sympy as sp
import numpy as np
from scipy.optimize import fsolve

def void_geometry(R_fuel,R_equiaxed,R_columnar,density_equiaxed_ratio,density_columnar_ratio.density_TD):
    #the main problem is the determination of equiaxed and columnar radius

    R_void = sp.symbols('R_void')
    R_void = sp.solvers.solve(R_fuel**2*density_TD-(R_columnar**2-R_void**2)*(density_columnar_ratio*density_TD)-(R_equiaxed**2-R_columnar**2)*(
            density_equiaxed_ratio*density_TD)-(R_fuel**2-R_equiaxed**2)*(density_TD))[0]

    return R_void 

  # Funzione per calcolare A in base a x e Pu_concentration
def A_f(x, Pu_concentration):
    A_value = 0.01926 + 1.06 * 10**-6 * x + 2.63 * 10**-8 * Pu_concentration
    return A_value

# Funzione per calcolare B in base a Pu_concentration
def B_f(Pu_concentration):
    B_value = 2.39 * 10**-4 + 1.37 * 10**-13 * Pu_concentration
    return B_value

# Funzione per calcolare la porosità in base alla densità
def porosity_f(density_value):
    p = 1 - density_value  # valutazione della porosità in diverse zone
    return p

# Definizione dei simboli
T = sp.symbols('T')
A = sp.symbols('A')
B = sp.symbols('B')
D = sp.symbols('D')
E = sp.symbols('E')
p = sp.symbols('p')
k = sp.Function('k')(T)

# Equazione da integrare per ottenere la conducibilità termica  
k_function = sp.Eq(k,1.755+(((1 / (A + B * T) + (D / T**2) * sp.exp(-E / T) * (1 - p)**2.5))-1.755)*sp.exp(-Burnup/128.75)) #il problema in questa equazione è il burnup
Equation_to_integrate = k_function.rhs

# Integrazione simbolica
integral_k = sp.integrate(Equation_to_integrate, T)

# Dati iniziali
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

# Creazione della matrice per i valori integrati di k
n = 1
m = 3
matrix_of_integrals = sp.MutableDenseMatrix(n, m, [0] * n * m)

# Stampa dell'integrale simbolico di k(T)
print("Symbolic integral of k(T) is:")
sp.pprint(integral_k)

# Calcolo e sostituzione dei valori nelle celle della matrice
for i in range(m):
    matrix_of_integrals[0, i] = integral_k.subs({
        A: A_value,
        B: B_value,
        D: D_value,
        E: E_value,
        p: porosity_f(DENSITY[i])
    })

# Stampa della matrice
print(matrix_of_integrals)

q_first_avg = 38.7 #kW/m

#evaluate all the equation to find the size of R_void, R_columnar, R_equiaxed

#I have a temperature profile which is different among all the z-interval(bottom, centre, top)
#Evaluate the case in the centre of the fuel element --> so considering z fixed and r can vary

T3 = 1000 - 273 #temperature max at the outer fuel radius (è una stima, non l'ho calcolata direttamente dalla funzione)
T2 = 1600 #temperature at the equiaxed radius
T1 = 1800 #temperature at the columnar radius
R_void = sp.symbols('R_void')
R_equiaxed = sp.symbols('R_equiaxed')
R_columnar = sp.symbols('R_columnar')
R_fuel = 0.542 #mm

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