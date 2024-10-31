def void_geometry(R_fuel,R_equiaxed,R_columnar,density_equiaxed_ratio,density_columnar_ratio.density_TD):
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

    # Definition of symbols
    T = sp.symbols('T')
    A = sp.symbols('A')
    B = sp.symbols('B')
    D = sp.symbols('D')
    E = sp.symbols('E')
    p = sp.symbols('p')
    k = sp.Function('k')(T)

    # Equation to integrate to obtain thermal conductivity
    k_function = sp.Eq(k, 1 / (A + B * T) + (D / T**2) * sp.exp(-E / T) * (1 - p)**2.5)
    Equation_to_integrate = k_function.rhs

    # Symbolic integration
    integral_k = sp.integrate(Equation_to_integrate, T)

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

    # Matrix creation for integrated k values
    n = 1
    m = 3
    matrix_of_integrals = sp.MutableDenseMatrix(n, m, [0] * n * m) #matrix contains all the integral_k equation for each zones

    # Print symbolic integral of k(T)
    print("Symbolic integral of k(T) is:")
    sp.pprint(integral_k)

    # Calculation and substitution of values into matrix cells
    
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