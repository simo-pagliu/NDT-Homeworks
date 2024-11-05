## attempt to build a function for Void Formation 

import numpy as np
import scipy.optimize as opt
import sympy as sp
import functions as f

def radius_finder(guess, beta, R_init, linear_heat_rate, Fuel_Properties, vars):

    h = 425 #mm
    T_R_toget = guess
    T_f_outer = f.get_temperature_at_point(h,R_init,vars.T_map)

    thermal_resistance = f.thermal_resistance_fuel(beta, Fuel_Properties)(T_f_outer)

    deltaT = linear_heat_rate/thermal_resistance

    R_toget = np.sqrt(1-(T_R_toget - T_f_outer)/deltaT)*R_init

    return R_toget, T_R_toget 

def R_void(Geometrical_Data, Fuel_Properties, R_col, R_eq):
    
    # from mass conservation equation using sympy
    R_void = sp.symbols('R_void')
    R_f_outer = sp.symbols('R_f_outer')
    R_columnar = sp.symbols('R_columnar')
    R_equiaxed = sp.symbols('R_equiaxed')
    density_af = sp.symbols('density_af')
    density_columnar = sp.symbols('density_columnar')
    density_equiaxed = sp.symbols('density_equiaxed')

    #write the equation LHS 
    LHS = R_f_outer**2*density_af
    RHS = (R_columnar**2-R_void**2)*density_columnar + (R_equiaxed**2-R_columnar**2)*density_equiaxed + (R_f_outer**2-R_equiaxed**2)*density_af

    # substitute the values of the variables
    LHS_sub = LHS.subs({R_f_outer: Geometrical_Data.fuel_outer_diameter/2, density_af: Fuel_Properties.Density})
    RHS_sub = RHS.subs({R_columnar: R_col, R_equiaxed: R_eq, density_columnar: Fuel_Properties.Theoretical_Density*0.98, density_equiaxed: Fuel_Properties.Theoretical_Density*0.95, R_f_outer: Geometrical_Data.fuel_outer_diameter/2, density_af: Fuel_Properties.Density})

    #solve the equation for R_void
    eq = sp.Eq(LHS_sub, RHS_sub)
    R_void = sp.solve(eq, R_void)[1]

    return R_void

    