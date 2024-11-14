##################################################
# Imports
##################################################
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import math
import dill
from IPython.display import display, Math
import nuclei_func as nf
import subprocess
import copy
##################################################
# Classes
##################################################
class Material_Proprieties:
    def __init__(self, Elements='', Qualities='', Density='',Theoretical_Density='',Percent_of_Theoretical_Density='', Emissivity ='', Molar_Mass='', Micro_Fission='', Micro_Absorption='', Viscosity='', Thermal_Conductivity='', Specific_Heat='', Thermal_Expansion_Coeff='', Melting_Temperature='', Boiling_Temperature='',Oxigen_to_metal_ratio='', Youngs_Modulus='', Poissons_Ratio='', Yield_Stress='', Ultimate_Tensile_Strength='', Nusselt_Number='', Grain_diameter=''):
        self.Elements = Elements
        self.Qualities = Qualities
        self.Theoretical_Density = Theoretical_Density
        self.Percent_of_Theoretical_Density = Percent_of_Theoretical_Density
        if Density == '':
            self.Density = self.Theoretical_Density * self.Percent_of_Theoretical_Density / 100
            self.Porosity = (self.Theoretical_Density - self.Density) / self.Density
        else:
            self.Density = Density
        self.Molar_Mass = Molar_Mass
        self.Micro_Fission = Micro_Fission
        self.Micro_Absorption = Micro_Absorption
        self.Viscosity = Viscosity
        self.Thermal_Conductivity = Thermal_Conductivity
        self.Specific_Heat = Specific_Heat
        self.Thermal_Expansion_Coeff = Thermal_Expansion_Coeff
        self.Emissivity = Emissivity
        self.Melting_Temperature = Melting_Temperature
        self.Boiling_Temperature = Boiling_Temperature
        self.Oxigen_to_metal_ratio = Oxigen_to_metal_ratio
        self.Youngs_Modulus = Youngs_Modulus
        self.Poissons_Ratio = Poissons_Ratio
        self.Yield_Stress = Yield_Stress
        self.Ultimate_Tensile_Strength = Ultimate_Tensile_Strength
        self.Nusselt_Number = Nusselt_Number
        self.Grain_diameter = Grain_diameter

        
# Geometry Data
class GeometryData:
    def __init__(self, fuel_outer_diameter, fuel_inner_diameter, cladding_outer_diameter, thickness_cladding, pin_pitch, h_values, fuel_pellet_height, fuel_roughness, cladding_roughness):
        self.fuel_outer_diameter = fuel_outer_diameter  # m
        self.fuel_inner_diameter = fuel_inner_diameter  # m (0 if solid fuel pellet)
        self.cladding_outer_diameter = cladding_outer_diameter  # m
        self.thickness_cladding = thickness_cladding  # m
        self.pin_pitch = pin_pitch # m
        self.h_values = h_values  # m
        self.fuel_pellet_height = fuel_pellet_height  # m
        self.nominal_size = [(cladding_outer_diameter - fuel_outer_diameter)/2 - thickness_cladding for cladding_outer_diameter, fuel_outer_diameter, thickness_cladding in zip(cladding_outer_diameter, fuel_outer_diameter, thickness_cladding)]
        self.effective_gap_size = [nominal_size + fuel_roughness + cladding_roughness for nominal_size in self.nominal_size]


# Thermo Hydraulics Specs
class ThermoHydraulicSpecs:
    def __init__(self, coolant_inlet_temp, coolant_inlet_pressure, coolant_mass_flow_rate, q_linear_avg, uptime, h_peak_factor, peak_factors, neutron_flux_peak):
        self.coolant_inlet_temp = coolant_inlet_temp  # K
        self.coolant_inlet_pressure = coolant_inlet_pressure  # Pa
        self.coolant_mass_flow_rate = coolant_mass_flow_rate  # kg/s
        self.q_linear_avg = q_linear_avg  # W/m
        self.uptime = uptime
        self.h_peak_factor = h_peak_factor  # W/m
        self.peak_factors = peak_factors
        self.neutron_flux_peak = neutron_flux_peak  # kg/s
        

# Define temperature map
class Temperature_Map:
    def __init__(self, X, Y, Z):
        self.r = X
        self.h = Y
        self.T = Z
        
# Define the class for Dimensioning Data
class DimensioningData:
    def __init__(self, filling_gas_pressure, filling_gas_temperature, temperature_map):
        self.p_gas = filling_gas_pressure  # Pa
        self.T_gas = filling_gas_temperature  # °C
        self.T_map = temperature_map

##################################################
# Radial values
##################################################
def create_meshgrid(geom_data):
    num_steps = len(geom_data.cladding_outer_diameter)
    r_values_list = []

    for i in range(num_steps):
        # Extract the corresponding values for this height
        cladding_outer_radius = geom_data.cladding_outer_diameter[i] / 2
        thickness_cladding = geom_data.thickness_cladding[i]
        fuel_outer_radius = geom_data.fuel_outer_diameter[i] / 2
        fuel_inner_radius = geom_data.fuel_inner_diameter[i] / 2

        # Calculate r_values for this specific height
        r_coolant_infinity = 10e-3 # 10 mm
        r_cladding_gap = cladding_outer_radius - thickness_cladding

        r_coolant = np.array([r_coolant_infinity])
        r_cladding = np.array([cladding_outer_radius])
        r_gap = np.linspace(r_cladding_gap, fuel_outer_radius, 4, endpoint=False)
        if fuel_inner_radius == 0:
            r_fuel = np.linspace(fuel_outer_radius, fuel_inner_radius, 25, endpoint=False)
            r_void = []
        else:
            r_fuel = np.linspace(fuel_outer_radius, fuel_inner_radius, 23, endpoint=False)
            r_void = [fuel_inner_radius, 0]

        # Concatenate all r_values for this height
        r_values_at_height = np.concatenate((r_coolant, r_cladding, r_gap, r_fuel, r_void))
        r_values_list.append(r_values_at_height)

    # Convert to a numpy array (each row corresponds to r_values at a specific height)
    r_values_array = np.array(r_values_list)
    h_values = geom_data.h_values

    # Create the meshgrid: X represents radius, Y represents height
    Y, X = np.meshgrid(h_values, np.arange(r_values_array.shape[1]), indexing='ij')

    # Replace X values with the corresponding r_values from r_values_array
    X = r_values_array

    return X, Y
##################################################
# Hydraulic Flow
##################################################
# Function to calculate the hydraulic flow parameters
def hydraulic_flow(thermo_hyd_spec, geom_data, coolant, temperature, i_h):
    passage_area = (1/2 * (geom_data.pin_pitch ** 2) * np.sin(np.pi/3)) - (1/2 * np.pi * (geom_data.cladding_outer_diameter[i_h]/2)**2)

    # Calculate the velocity of the fluid
    velocity = thermo_hyd_spec.coolant_mass_flow_rate / (coolant.Density(temperature) * passage_area)

    # Calculate the hydraulic diameter
    wetted_perimeter = np.pi * geom_data.cladding_outer_diameter[i_h] / 2

    hydraulic_diameter = (4 * passage_area) / (wetted_perimeter)

    return velocity, passage_area, hydraulic_diameter

##################################################
# Coolant Thermohydraulics
##################################################
def coolant_thermohydraulics(geom_data, thermo_hyd_spec, coolant, temperature, i_h):
    density = coolant.Density(temperature)
    viscosity = coolant.Viscosity(temperature)
    thermal_conductivity = coolant.Thermal_Conductivity(temperature)
    c_p = coolant.Specific_Heat(temperature)

    # Calculate the velocity and passage area
    velocity, _, d_h = hydraulic_flow(thermo_hyd_spec, geom_data, coolant, temperature, i_h)

    # Adimensional numbers
    reynolds = (density * velocity * d_h) / viscosity
    prandtl = c_p * viscosity / thermal_conductivity
    peclet = reynolds * prandtl
    nusselt = coolant.Nusselt_Number(peclet)

    # HTC calculation
    htc = nusselt * thermal_conductivity / d_h

    # Calculate the thermal resistance
    radius_cladding_out = geom_data.cladding_outer_diameter[i_h] / 2
    thermal_resistance = 1 / (2 * np.pi * radius_cladding_out * htc)

    return thermal_resistance, velocity

##################################################
# Power and Neutron FLux Profile
##################################################
def power_profile(h, thermo_hyd_spec, value = 'power'):
    pf = thermo_hyd_spec.peak_factors
    h_from_c = thermo_hyd_spec.h_peak_factor
    if value == 'power':
        q = thermo_hyd_spec.q_linear_avg
        peak_value = q * len(pf) / sum(pf)
    elif value == 'neutron_flux':
        peak_value = thermo_hyd_spec.neutron_flux_peak
    # Computes peak power such that the average power is q_linear_avg

    # Compute power values for each interval
    power_values = [peak_value * factor for factor in pf]

    # Find the interval that h belongs to
    length_h_from_c = len(h_from_c)
    interval_boundaries = [0] + [(h_from_c[i] + h_from_c[i + 1]) / 2 for i in range(length_h_from_c - 1)] + [850]

    for i in range(length_h_from_c):
        if interval_boundaries[i] <= h <= interval_boundaries[i + 1]:
            return power_values[i]

##################################################
# Thermal Resistances
##################################################
def thermal_resistance_cladding(geom_data, cladding, temperature, i_h):
    radius_cladding_out = geom_data.cladding_outer_diameter[i_h] / 2
    radius_cladding_in = radius_cladding_out - geom_data.thickness_cladding[i_h]
    k = cladding.Thermal_Conductivity(temperature)
    thermal_resistance = np.log(radius_cladding_out / radius_cladding_in) / (2 * np.pi * k)
    return thermal_resistance

def thermal_resistance_gap(geom_data, helium, fuel, cladding, t_gas, He_percentage, t_fuel_out, i_h):
    radius_gap_in = geom_data.fuel_outer_diameter[i_h] / 2

    # Conductive heat transfer
    conduction = helium.Thermal_Conductivity(t_gas, He_percentage) / geom_data.effective_gap_size[i_h]

    # Radiative heat transfer
    radiation = 4 * 5.67e-8 * (t_fuel_out**3) / (1/(fuel.Emissivity) + 1/(cladding.Emissivity) - 1)

    # Overall heat transfer coefficient
    htc = conduction + radiation

    # Calculate the thermal resistance of the gap
    thermal_resistance = 1 / (2 * np.pi * radius_gap_in * htc)

    return thermal_resistance

def thermal_resistance_fuel(Burnup, fuel, temperature):
    A = 0.01926 + 1.06e-6 * fuel.Oxigen_to_metal_ratio + 2.63e-8 * fuel.Molar_Mass[-1]
    B = 2.39e-4 + 1.37e-13 * fuel.Molar_Mass[-1]
    D = 5.27e9
    E = 17109.5
    
    # Accounts for fuel regions due to restructuring
    if temperature > 1800:
        Porosity = 0.00
    elif temperature > 1600:
        Porosity = 0.02
    else:
        Porosity = fuel.Porosity
        
    # Calculate k_0
    k_0 = (1 / (A + B * temperature) + D / (temperature**2) * np.exp(-E / temperature)) * (1 - Porosity)**2.5

    # Calculate final thermal conductivity k
    k = 1.755 + (k_0 - 1.755) * np.exp(-Burnup / 128.75)

    # Calculate the thermal resistance of the fuel
    thermal_resistance = 1 / (4 * np.pi * k)
    return thermal_resistance

##################################################
# Temperature Map
##################################################
def temperature_map(coolant, cladding, gap, fuel, thermo_hyd_spec, geom_data, T_fuel_out, Burnup, He_percentage):
    h_values = geom_data.h_values
    # r_values = calculate_r_values(geom_data)

    # Compute the height of each slice
    dz = h_values[1] - h_values[0]
    
    # Create the meshgrid
    X, Y = create_meshgrid(geom_data)
    Z = []
    
    # Initialize the temperature
    Temp_coolant = thermo_hyd_spec.coolant_inlet_temp

    # Compute the temperature profile for each height step
    for i_h, h in enumerate(h_values): 
        # Get radii
        r_values = X[i_h, :]

        # Compute temperature of the coolant
        f = 1/2
        c_p = coolant.Specific_Heat(Temp_coolant)
        m_rate = thermo_hyd_spec.coolant_mass_flow_rate

        power = power_profile(h, thermo_hyd_spec)

        Temp_coolant = Temp_coolant + power * dz / ( m_rate / f * c_p)
        
        # Initialize the temperature profile
        T_radial = [Temp_coolant] # Temperature of the coolant (ideally at r = infinity)

        # Create Limits
        r_coolant_cladding = geom_data.cladding_outer_diameter[i_h] / 2
        r_cladding_gap = geom_data.cladding_outer_diameter[i_h] / 2 - geom_data.thickness_cladding[i_h]
        r_gap_fuel = geom_data.fuel_outer_diameter[i_h] / 2
        r_fuel_in = geom_data.fuel_inner_diameter[i_h] / 2

        # Index of the interfaces
        idx_fuel = np.argmin(np.abs(r_values - r_gap_fuel))
        idx_gap = np.argmin(np.abs(r_values - r_cladding_gap))

        # Compute the temperature profile
        for j, r in enumerate(r_values[1:], start=1): 
            dr = r_values[j-1] - r_values[j]

            # In the void
            if r < r_fuel_in:
                T_value = T_radial[j-1]
    
            # In the fuel
            elif r < r_gap_fuel:
                
                # Colunmar Region
                if T_radial[j-1] > 1800:
                    fuel.Percent_of_Theoretical_Density = 100
                    
                # Equiaxed Region
                elif T_radial[j-1] > 1600:
                    fuel.Percent_of_Theoretical_Density = 98

                # Thermal resistance of the fuel
                th_res = thermal_resistance_fuel(Burnup, fuel, T_radial[j-1]) 
                # Compute the temperature
                T_value = T_radial[idx_fuel] + power * th_res * (1 - (r/r_gap_fuel)**2)
            
            # In the gap
            elif r < r_cladding_gap:
                # Thermal resistance of the gap
                th_res = thermal_resistance_gap(geom_data, gap, fuel, cladding, T_radial[j-1], He_percentage, T_fuel_out, i_h)
                # Compute the temperature
                T_value = T_radial[idx_gap] + power * th_res * np.log(r_cladding_gap / r)
                
            # In the cladding
            elif r < r_coolant_cladding:
                # Thermal resistance of the cladding
                th_res = thermal_resistance_cladding(geom_data, cladding, T_radial[j-1], i_h)
                # Compute the temperature
                T_value = T_radial[j-1] + power * th_res * (dr / (r_coolant_cladding - r_cladding_gap))
                
            # In the coolant
            else:
                # Thermal resistance of the coolant
                th_res, _ = coolant_thermohydraulics(geom_data, thermo_hyd_spec, coolant, T_radial[0], i_h)
                # Compute the temperature
                T_value = T_radial[j-1] + power * th_res * (dr / (r_values[0] - r_coolant_cladding))

            # Compute new value of T
            T_radial.append(T_value)

        Z.append(T_radial)

    Z = np.array(Z)
    T_map = Temperature_Map(X, Y, Z)
    return T_map

##################################################
# Search Functions
##################################################
def get_temperature_at_point(h_requested, r_requested,T_map):
    h_values = T_map.h
    r_values = T_map.r
    T_values = T_map.T
    h_idx = np.argmin(np.abs(h_values[:, 0] - h_requested))
    r_idx = np.argmin(np.abs(r_values[0, :] - r_requested))
    return T_values[h_idx, r_idx]

def get_radius_at_temperature(T_requested,T_map):
    r_values = []
    T_values = T_map.T
    for t,r in zip(T_values,T_map.r):
        T_idx = np.argmin(np.abs(t - T_requested))
        r_values.append(r[T_idx])
    
    return r_values

##################################################
# Void Swelling
##################################################
def void_swelling(T_map, geom_data, thermo_hyd_spec):
    Volume_expansion_fission_gas = []
    
    for idx_h, h in enumerate(T_map.h[:, 0]):
        r_fuel_out = geom_data.fuel_outer_diameter[idx_h] / 2
        idx_fuel = np.argmin(np.abs(T_map.r[idx_h, :] - r_fuel_out))
        r_vals = T_map.r[idx_h, idx_fuel:-1]

        phi = power_profile(h, thermo_hyd_spec, value = 'neutron_flux') * thermo_hyd_spec.uptime

        temperature = [get_temperature_at_point(h, r, T_map) for r in r_vals]
        temperature_avg = np.mean(temperature)

        temp = 1.5e-3 * np.exp(-2.5 * ((temperature_avg - 273 - 450) / 100) ** 2) * (phi / 1e22) ** 2.75
        Volume_expansion_fission_gas.append(temp)
    return Volume_expansion_fission_gas

def get_R_void(Fuel_Proprieties, R_col, R_eq):
        
    R_void = []
    density_af = Fuel_Proprieties.Density
    density_columnar=Fuel_Proprieties.Theoretical_Density*0.98
    density_equiaxed= Fuel_Proprieties.Theoretical_Density*0.95
    
    for r in range(len(R_col)):
        
        R_voidd = np.sqrt(R_col[r]**2 - R_eq[r]**2*(density_af/density_columnar) +(R_eq[r]**2 - R_col[r]**2)*(density_equiaxed/density_af))
        R_void.append(R_voidd)
        
    return R_void

##################################################
# Thermal Expansion
##################################################
def thermal_expansion(fuel, cladding, cold_geometrical_data, T_map):
    # Initialize
    h_vals = cold_geometrical_data.h_values
    
    #### Fuel ####
    # Initialize
    fuel_inner_hot = []
    fuel_outer_hot = []
    # Reference Temperatrue
    T_0 = 25 + 273.15
    # Calculate for fuel
    fuel_alpha = fuel.Thermal_Expansion_Coeff

    # Compute along the height
    for i_h, h in enumerate(h_vals):
        # Starting values
        fuel_outer = cold_geometrical_data.fuel_outer_diameter[i_h] / 2
        fuel_inner = cold_geometrical_data.fuel_inner_diameter[i_h] / 2
        # Compute along the radius
        for r in [fuel_inner, fuel_outer]:
            T_now = get_temperature_at_point(h, r, T_map)
            Expanded_Radius = r * (1 + fuel_alpha * (T_now - T_0))
            if r == fuel_inner:
                fuel_inner_hot.append(2*Expanded_Radius)
            else:
                fuel_outer_hot.append(2*Expanded_Radius)

    #### Cladding ####
    # Initialize
    thickness_cladding = []
    cladding_outer_hot = []
    # Expansion coefficient
    cladding_alpha = cladding.Thermal_Expansion_Coeff

    # Compute along the height
    for i_h, h in enumerate(h_vals):
        # Starting values
        cladding_outer_radius = cold_geometrical_data.cladding_outer_diameter[i_h] / 2
        cladding_inner_radius = cladding_outer_radius - cold_geometrical_data.thickness_cladding[i_h]
        # Compute along the radius
        for r in [cladding_outer_radius, cladding_inner_radius]:
            T_now = get_temperature_at_point(h, r, T_map)
            strain = cladding_alpha(T_now)
            Expanded_Radius = r * (1 + strain)
            if r == cladding_outer_radius:
                # Outer radius computation
                cladding_outer_hot.append(2*Expanded_Radius)
            else:
                # Inner radius computation
                thickness = cladding_outer_hot[i_h]/2 - Expanded_Radius
                thickness_cladding.append(thickness)

    return fuel_outer_hot, fuel_inner_hot, cladding_outer_hot, thickness_cladding

##################################################
# Fission Gas Production
##################################################
# Size of bubble radius
a = 10e-6  # m

# Starting gas temperature
T_gas = 20 + 273 # °C --> K

# Initial gas pressure
p_gas = 1e5  # Pa

# Fission cross sections (see file "Useful Data.xlsx")
sigma_235 = 1.047756375  # barn
sigma_238 = 0.55801001  # barn
sigma_pu = 1.689844625  # barn

# Fission yield
fission_yield = 0.3

# Diffusivity evaluation parameters [Matzke, 1980]
d_0 = 5e-8  # m^2/s
q = 40262

def fission_gas_production(h_plenum, Fuel_Proprieties, ThermoHydraulics, Geometrical_Data, T_map):
    """
    Function to compute the rate theory fission gas calculations and generate relevant plots.
    Outputs:
        He_percentage: Percentage of helium present inside the fuel still trapped in
        new_p_gas: Gas pressure in the plenum
    """
    # Calculate molar mass
    molar_mass = nf.mixture(Fuel_Proprieties.Molar_Mass, Fuel_Proprieties.Qualities, normalization_cond='normalize')  # g/mol

    # Calculate macroscopic cross sections
    macro_235 = nf.macro(sigma_235, Fuel_Proprieties.Density, molar_mass)  # cm^-1
    macro_238 = nf.macro(sigma_238, Fuel_Proprieties.Density, molar_mass)  # cm^-1
    macro_pu = nf.macro(sigma_pu, Fuel_Proprieties.Density, molar_mass)  # cm^-1

    fission_xs = nf.mixture([macro_235, macro_238, 0, macro_pu], Fuel_Proprieties.Qualities)  # cm^-1

    # Compute average neutron flux
    # power profile
    q_values = [power_profile(h, ThermoHydraulics) for h in Geometrical_Data.h_values]
    peak_to_average_ratio = max(q_values) / ThermoHydraulics.q_linear_avg
    average_neutron_flux = ThermoHydraulics.neutron_flux_peak / peak_to_average_ratio

    # Calculate average fission rate
    avg_fission_rate = average_neutron_flux * fission_xs * 1e6  # [fissions/m^3 s]

    # Diffusivity function
    diffusivity = lambda temperature: d_0 * np.exp(-q / temperature)  # m^2/s

    # Compute the average maximum fuel temperature
    points = [0.425, 0, 0.850]
    fuel_inner_diameter_avg = np.mean(Geometrical_Data.fuel_inner_diameter)
    temperature_max_average_fuel = sum(
        get_temperature_at_point(point, fuel_inner_diameter_avg / 2, T_map)
        for point in points
    ) / len(points)

    diffusivity_coeff = diffusivity(temperature_max_average_fuel)  # m^2/s

    # Define P as a lambda function based on the given solution
    P_lambda = lambda t: fission_yield * avg_fission_rate * t

    # Define GM as a lambda function based on the given solution
    GM_lambda = lambda r: (-fission_yield * avg_fission_rate * r**2 / (4 * diffusivity_coeff)) + (2.5e-11 * fission_yield * avg_fission_rate / diffusivity_coeff)

    ### Integrate the solution for GM in the radial direction using numerical method
    r_vals = np.linspace(0, a, 1000)  # Create a range of radius values
    GM_vals = [GM_lambda(r) for r in r_vals]  # Evaluate GM_lambda at each radius value
    GM_final = np.trapz(GM_vals, r_vals)  # Numerically integrate using trapezoidal rule

    ## Compute the amount of gas released in plenum after 1 year

    # Production of fission gas in the fuel
    time = 360 * 24 * 3600  # 1 year (of operation) in seconds
    total_fission_gas = P_lambda(time)  # Total amount of fission gas produced

    # Total amount of fission gas inside the grains
    n_grains_pellet = 1e5  # Number of grains in a pellet
    n_pellets_pin = round(850e-3 / Geometrical_Data.fuel_pellet_height)  # Number of pellets in a pin
    total_fission_gas_grains = GM_final * n_grains_pellet * n_pellets_pin

    # Total amount of fission gas released in the plenum
    total_fission_gas_released = total_fission_gas - total_fission_gas_grains

    # Calculate the percentage of helium trapped inside the fuel
    He_percentage = (total_fission_gas_grains / total_fission_gas)

    # Vector of possible plenum heights
    fuel_column_height = 850e-3  # m

    # Corresponding volume to accommodate gases
    r_cladding_gap = np.mean(np.array(Geometrical_Data.cladding_outer_diameter) / 2 - np.array(Geometrical_Data.thickness_cladding))
    r_gap_fuel = np.mean(np.array(Geometrical_Data.fuel_outer_diameter) / 2)
    
    V_plenum = (np.pi * r_cladding_gap**2 * h_plenum) + (np.pi * (r_cladding_gap**2 - r_gap_fuel**2) * fuel_column_height)

    # Find initial quantity of He present in the plenum
    initial_moles_he = p_gas * V_plenum / (8.314 * (T_gas))  # moles

    # Find the additional moles of fission gases released in the plenum
    fuel_outer_diameter_avg = np.mean(Geometrical_Data.fuel_outer_diameter)
    V_pin = np.pi * (fuel_outer_diameter_avg / 2)**2 * fuel_column_height  # m^3
    additional_moles_fg = (total_fission_gas_released * V_pin) / 6.022e23  # moles

    # Find the total moles of gases in the plenum
    total_moles_gas = initial_moles_he + additional_moles_fg  # moles

    # Find the new pressure in the plenum
    new_p_gas = total_moles_gas * 8.314 * (T_gas) / V_plenum  # Pa

    return He_percentage, new_p_gas

##################################################
# Main Function
##################################################
def initialize_params():
    """
    Initialization function to set up parameters.
    
    Returns:
        dict: Dictionary containing initialized parameters.
    """
    # Load all objects in a single line using dill
    with open("data_file.dill", "rb") as file:
        data_objects = dill.load(file)

    params = {
        "Cladding_Proprieties": data_objects["Cladding_Proprieties"],
        "Fuel_Proprieties": data_objects["Fuel_Proprieties"],
        "Coolant_Proprieties": data_objects["Coolant_Proprieties"],
        "Helium_Proprieties": data_objects["Helium_Proprieties"],
        "ThermoHydraulics": data_objects["ThermoHydraulics"],
        "Geometrical_Data_Cold": data_objects["Geometrical_Data"],
    }
    return params


def compute_burnup(params):
    """
    Compute the burnup value.
    
    Args:
        params (dict): Dictionary containing parameters.
    
    Returns:
        float: Computed burnup value.
    """
    Geometrical_Data_Cold = params["Geometrical_Data_Cold"]
    Fuel_Proprieties = params["Fuel_Proprieties"]
    ThermoHydraulics = params["ThermoHydraulics"]

    MOX_volume = np.pi * (Geometrical_Data_Cold.fuel_outer_diameter[0]/2)**2 * Geometrical_Data_Cold.h_values[-1]
    MOX_mass = MOX_volume * Fuel_Proprieties.Density * 1000
    MOX_molar_mass = nf.mixture(Fuel_Proprieties.Molar_Mass, Fuel_Proprieties.Qualities) * 1e-3  # g/mol --> kg/mol
    MOX_mols = MOX_mass / MOX_molar_mass
    HM_molar_mass = [Fuel_Proprieties.Molar_Mass[i] for i in [0, 1, 3]]
    HM_qualities = [Fuel_Proprieties.Qualities[i] for i in [0, 1, 3]]
    HM_molar_mass = nf.mixture(HM_molar_mass, HM_qualities) * 1e-3 # g/mol --> kg/mol
    HM_mass = MOX_mols * HM_molar_mass  # kg of Heavy Metal

    Burnup = ThermoHydraulics.uptime / (24 * 3600) \
        * ThermoHydraulics.q_linear_avg * Geometrical_Data_Cold.h_values[-1] \
        / HM_mass
    Burnup = Burnup * 1e-6  # Wd/kgU * 1e-9/1e-3 --> GWd/t(HM)
    return Burnup


def reset_geometrical_data(params, delta):
    """
    Reset geometrical data with a new cladding thickness.
    
    Args:
        params (dict): Dictionary containing parameters.
        delta (float): New cladding thickness.
    
    Returns:
        Geometrical_Data: Reset Geometrical_Data object.
    """
    Geometrical_Data = dill.loads(dill.dumps(params["Geometrical_Data_Cold"]))  # Reset Geometrical_Data
    Geometrical_Data.thickness_cladding = np.full_like(Geometrical_Data.h_values, delta)
    return Geometrical_Data


def update_temperatures(params, Geometrical_Data, T_fuel_out, Burnup, He_percentage, h_plenum, previous_T_map):
    """
    Update temperature map and geometrical data.
    
    Args:
        params (dict): Dictionary containing parameters.
        Geometrical_Data: Geometrical_Data object.
        T_fuel_out (float): Initial guess for fuel outer temperature.
        Burnup (float): Burnup value.
    
    Returns:
        tuple: Updated temperature map, fuel outer temperature, and Geometrical_Data.
    """

    # Thermal Expansion
    Geometrical_Data.fuel_outer_diameter, \
    Geometrical_Data.fuel_inner_diameter, \
    Geometrical_Data.cladding_outer_diameter, \
    Geometrical_Data.thickness_cladding = thermal_expansion(params["Fuel_Proprieties"], params["Cladding_Proprieties"], params["Geometrical_Data_Cold"], previous_T_map)
    
    # Void Formation due to Restructuring
    R_equiaxied = get_radius_at_temperature(1600, previous_T_map)
    R_columnar = get_radius_at_temperature(1800, previous_T_map)
    R_void = get_R_void(params["Fuel_Proprieties"], R_columnar, R_equiaxied)
    Geometrical_Data.fuel_inner_diameter = [2 * r for r in R_void]

    # Void Swelling
    void_swell = void_swelling(previous_T_map, Geometrical_Data, params["ThermoHydraulics"])
    radius_before = [d/ 2 for d in Geometrical_Data.fuel_outer_diameter]
    radius_swelled = [radius_before * np.sqrt(1 + swell) for radius_before, swell in zip(radius_before, void_swell)]
    Geometrical_Data.fuel_outer_diameter = [2 * r for r in radius_swelled]

    # Fission Gas Production
    He_percentage, Gas_Pressure = fission_gas_production(h_plenum, params["Fuel_Proprieties"], params["ThermoHydraulics"], Geometrical_Data, previous_T_map)

    T_map = temperature_map(params["Coolant_Proprieties"], params["Cladding_Proprieties"], params["Helium_Proprieties"], 
                              params["Fuel_Proprieties"], params["ThermoHydraulics"], Geometrical_Data, T_fuel_out, Burnup, He_percentage)
    idx_fuel = np.argmin(np.abs(T_map.r[5, :] - Geometrical_Data.fuel_outer_diameter[0]/2))
    T_fuel_out = T_map.T[5, idx_fuel]

    return T_map, T_fuel_out, Geometrical_Data, He_percentage, Gas_Pressure


def main(params):
    """
    Main function that takes initialized parameters and executes the main logic.
    
    Args:
        params (dict): Dictionary containing parameters.
    """
    Burnup = compute_burnup(params)
    print("Burnup:", Burnup)

    # VARIABLES TO OPTIMIZE
    cladding_thicknesses = [50e-6] #np.linspace(565e-6, 1e-6, 6)
    h_plenum = 1.8 # m

    # Initial values
    T_fuel_out = 1160  # K (Initial guess)
    He_percentage = 1  # Initial Value
    iterations_gap = []
    iterations_plenum = []

    # Loop over cladding thicknesses
    for delta in cladding_thicknesses:
        residual = 1  # Placeholder for residual
        j = 0
        Geometrical_Data = reset_geometrical_data(params, delta)
        previous_T_map = temperature_map(params["Coolant_Proprieties"], params["Cladding_Proprieties"], params["Helium_Proprieties"], 
                              params["Fuel_Proprieties"], params["ThermoHydraulics"], Geometrical_Data, T_fuel_out, Burnup, He_percentage)

        while residual > 1e-3:
            j += 1

            r_cladding_gap = np.array(Geometrical_Data.cladding_outer_diameter) / 2 - np.array(Geometrical_Data.thickness_cladding)
            r_gap_fuel = np.array(Geometrical_Data.fuel_outer_diameter) / 2
            diff = r_cladding_gap - r_gap_fuel

            if np.all(diff < 0):
                print(f"\033[91mFuel Diameter < Cladding Diameter, Gap is closed at iteration {j}\033[0m")
                break
            else:
                T_map, T_fuel_out, Geometrical_Data, He_percentage, Gas_Pressure = update_temperatures(params, Geometrical_Data, T_fuel_out, Burnup, He_percentage, h_plenum, previous_T_map)

                residual = np.mean(np.abs(T_map.T - previous_T_map.T)) / np.mean(previous_T_map.T)
                previous_T_map = copy.deepcopy(T_map)

        iterations_gap.append(T_map)
        iterations_plenum.append(Gas_Pressure)

    plt.figure()
    for T_map, delta, Gas_p in zip(iterations_gap, cladding_thicknesses, iterations_plenum):
        plt.plot(T_map.r[0, :], T_map.T[0, :], label=f"{delta*1e3:.2f} mm ---> Gap at {Gas_p*1e-6:.2f} MPa", marker='o')
    plt.axvline(x=Geometrical_Data.fuel_outer_diameter[0]/2, color='r', linestyle='--', label='Fuel Outer Diameter')
    plt.axvline(x=Geometrical_Data.fuel_inner_diameter[0]/2, color='g', linestyle='--', label='Fuel Inner Diameter')
    plt.axvline(x=Geometrical_Data.cladding_outer_diameter[0]/2, color='b', linestyle='--', label='Cladding Outer Diameter')
    plt.axvline(x=Geometrical_Data.cladding_outer_diameter[0]/2 - Geometrical_Data.thickness_cladding[0], color='y', linestyle='--', label='Cladding Inner Diameter')
    plt.title("Temperature Profile")
    plt.xlabel("Radius [m]")
    plt.ylabel("Temperature [K]")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    params = initialize_params()
    print("Parameters initialized.")
    main(params)