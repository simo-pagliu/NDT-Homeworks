## attempt to build a function for Void Formation 

import numpy as np
import scipy.optimize as opt
import sympy as sp
import functions as f

# def radius_finder(guess, beta, R_init, linear_heat_rate, Fuel_Properties, vars):

#     h = 425 #mm
#     T_R_toget = guess
#     T_f_outer = f.get_temperature_at_point(h,R_init,vars.T_map)

#     thermal_resistance = f.thermal_resistance_fuel(beta, Fuel_Properties)(T_f_outer)

#     deltaT = linear_heat_rate/thermal_resistance

#     R_toget = np.sqrt(1-(T_R_toget - T_f_outer)/deltaT)*R_init

#     return R_toget, T_R_toget 
##################################################
# Central Void Formation
##################################################
# def get_radius_at_temprature(T_requested,T_map):

#     T_values = T_map.T
#     T_idx = np.argmin(np.abs(T_values - T_requested))
#     print(T_idx)
#     return T_idx

# xxx = get_radius_at_temprature(1800, T_map)

#######
def get_radius_at_temprature(T_requested,T_map):

    r_values = []
    T_values = T_map.T
    for t,r,h in zip(T_values,T_map.r,T_map.h):
        T_idx = np.argmin(np.abs(t - T_requested))
        r_values.append(r[T_idx])
    
    return r_values

def get_R_void(Fuel_Properties, R_col, R_eq):
        
    R_void = []
    density_af = Fuel_Properties.Density
    density_columnar=Fuel_Properties.Theoretical_Density*0.98
    density_equiaxed= Fuel_Properties.Theoretical_Density*0.95
    
    for r in range(len(R_col)):
        
        R_voidd = np.sqrt(R_col[r]**2 - R_eq[r]**2*(density_af/density_columnar) +(R_eq[r]**2 - R_col[r]**2)*(density_equiaxed/density_af))
        R_void.append(R_voidd)
        
    return R_void

    