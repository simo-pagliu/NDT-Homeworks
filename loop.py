##################################################
# Imports
##################################################
import numpy as np
# import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# import plotly.graph_objects as go
import dill
import nuclei_func as nf
from scipy.interpolate import interp1d
# import subprocess
import math # Used in some data
import copy
##################################################
# Classes
##################################################
class Material_Proprieties:
    def __init__(self, Elements='', Qualities='', Density='',Theoretical_Density='',Percent_of_Theoretical_Density='', Porosity_Columnar='', Porosity_Equiaxed='', Emissivity ='', Molar_Mass='', Micro_Fission='', Micro_Absorption='', Viscosity='', Thermal_Conductivity='', Specific_Heat='', Thermal_Expansion_Coeff='', Melting_Temperature='', Boiling_Temperature='',Oxigen_to_metal_ratio='', Youngs_Modulus='', Poissons_Ratio='', Yield_Stress='', Ultimate_Tensile_Strength='', Nusselt_Number='', Grain_diameter='', Starting_Gas_Temperature='', Initial_Gas_Pressure='', Sigma_235='', Sigma_238='', Sigma_Pu='', Fission_Yield=''):
        self.Elements = Elements
        self.Qualities = Qualities
        self.Theoretical_Density = Theoretical_Density
        self.Percent_of_Theoretical_Density = Percent_of_Theoretical_Density
        if Density == '':
            self.Density = self.Theoretical_Density * self.Percent_of_Theoretical_Density / 100
            self.Porosity = (self.Theoretical_Density - self.Density) / self.Density
        else:
            self.Density = Density
        self.Porosity_Columnar = Porosity_Columnar
        self.Porosity_Equiaxed = Porosity_Equiaxed
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
        self.Starting_Gas_Temperature = Starting_Gas_Temperature
        self.Initial_Gas_Pressure = Initial_Gas_Pressure
        self.Sigma_235 = Sigma_235
        self.Sigma_238 = Sigma_238
        self.Sigma_Pu = Sigma_Pu
        self.Fission_Yield = Fission_Yield
        self.Nusselt_Number = Nusselt_Number
        self.Grain_diameter = Grain_diameter
        
# Geometry Data
class GeometryData:
    def __init__(self, fuel_outer_diameter, fuel_inner_diameter, cladding_outer_diameter, thickness_cladding, pin_pitch, h_values, fuel_pellet_height, fuel_roughness, cladding_roughness, Initial_Gas_Temperature='', Initial_Gas_Pressure='', Sigma_235='', Sigma_238='', Sigma_Pu='', Fission_Yield=''):
        self.Initial_Gas_Temperature = Initial_Gas_Temperature
        self.Initial_Gas_Pressure = Initial_Gas_Pressure
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
    def __init__(self, coolant_inlet_temp, coolant_inlet_pressure, coolant_mass_flow_rate, q_linear_avg, uptime, h_peak_factor, peak_factors, neutron_flux_peak, Starting_Gas_Temperature='', Initial_Gas_Pressure='', Sigma_235='', Sigma_238='', Sigma_Pu='', Fission_Yield=''):
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
    def __init__(self, X, Y, Z, Starting_Gas_Temperature='', Initial_Gas_Pressure='', Sigma_235='', Sigma_238='', Sigma_Pu='', Fission_Yield=''):
        self.r = X
        self.h = Y
        self.T = Z
        
# Define the class for Dimensioning Data
class DimensioningData:
    def __init__(self, filling_gas_pressure, filling_gas_temperature, temperature_map, Starting_Gas_Temperature='', Initial_Gas_Pressure='', Sigma_235='', Sigma_238='', Sigma_Pu='', Fission_Yield=''):
        self.p_gas = filling_gas_pressure  # Pa
        self.T_gas = filling_gas_temperature  # °C
        self.T_map = temperature_map

coolant_infinity_limit = 6e-3 # mm
##################################################
# Radial values
##################################################
def create_meshgrid(geom_data):
    r_values_list = []
    for i in range(len(geom_data.h_values)):
        # Extract the corresponding values for this height
        cladding_outer_radius = geom_data.cladding_outer_diameter[i] / 2
        thickness_cladding = geom_data.thickness_cladding[i]
        fuel_outer_radius = geom_data.fuel_outer_diameter[i] / 2
        fuel_inner_radius = geom_data.fuel_inner_diameter[i] / 2

        # Calculate r_values for this specific height
        r_coolant_infinity = cladding_outer_radius + coolant_infinity_limit  # 8 mm
        r_cladding_gap = cladding_outer_radius - thickness_cladding

        coolant_points = 2
        r_coolant = np.linspace(r_coolant_infinity, cladding_outer_radius, coolant_points, endpoint=False)

        gap_points = 5
        cladding_points = 2
        available_gap = cladding_outer_radius - thickness_cladding - fuel_outer_radius
        if available_gap > 0:
            r_gap = np.linspace(r_cladding_gap, fuel_outer_radius, gap_points, endpoint=False)
            r_cladding = np.linspace(cladding_outer_radius, r_cladding_gap, cladding_points, endpoint=False)
        else:
            r_gap = []
            r_cladding = np.linspace(cladding_outer_radius, fuel_outer_radius, cladding_points + gap_points, endpoint=False)


        fuel_points = 40
        if fuel_inner_radius == 0:
            r_fuel = np.linspace(fuel_outer_radius, fuel_inner_radius, fuel_points, endpoint=False)
            r_void = []
        else:
            r_fuel = np.linspace(fuel_outer_radius, fuel_inner_radius, fuel_points - 2, endpoint=False)
            r_void = [fuel_inner_radius, 0]

        # Concatenate all r_values for this height
        r_values_at_height = np.concatenate((r_coolant, r_cladding, r_gap, r_fuel, r_void))
        r_values_list.append(r_values_at_height)

    # Convert to a numpy array (each row corresponds to r_values at a specific height)
    r_values_array = np.array(r_values_list, dtype=object)
    h_values = geom_data.h_values

    # Create the meshgrid: X represents radius, Y represents height
    # max_len = max(len(row) for row in r_values_array)
    # padded_r_values_array = np.array([np.pad(row, (0, max_len - len(row)), constant_values=np.nan) for row in r_values_array])
    X, Y = np.meshgrid(r_values_array[0], h_values)
    X = r_values_array

    # Replace X values with the corresponding r_values from padded_r_values_array
    # X = padded_r_values_array

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
        # q = thermo_hyd_spec.q_linear_avg
        # peak_value = q * len(pf) / sum(pf)
        peak_value = thermo_hyd_spec.q_linear_avg
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
    gap_size = max(geom_data.effective_gap_size[i_h], 1e-6)  # Prevent zero or negative gap size
    conduction = helium.Thermal_Conductivity(t_gas, He_percentage) / gap_size

    # Radiative heat transfer
    radiation = 4 * 5.67e-8 * (t_fuel_out**3) / (1/(fuel.Emissivity) + 1/(cladding.Emissivity) - 1)

    # Overall heat transfer coefficient
    htc = conduction + radiation

    # Calculate the thermal resistance of the gap
    thermal_resistance = 1 / (2 * np.pi * radius_gap_in * htc)

    return thermal_resistance

def thermal_resistance_fuel(Burnup, fuel, temperature):
    pu_wt = nf.w2mol([fuel.Molar_Mass[-1]], [239])[0]
    A = 0.01926 + 1.06e-6 * fuel.Oxigen_to_metal_ratio + 2.63e-8 * pu_wt
    B = 2.39e-4 + 1.37e-13 * pu_wt
    D = 5.27e9
    E = 17109.5
    
    # Accounts for fuel regions due to restructuring
    if temperature > 1800:
        Porosity = fuel.Porosity_Columnar
    elif temperature > 1600:
        Porosity = fuel.Porosity_Equiaxed
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

    # Compute the height of each slice
    dz = h_values[1] - h_values[0]
    
    # Create the meshgrid
    X, Y = create_meshgrid(geom_data)
    Z = []
    velocity_vector = [] # Velocity vector
    
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

        Temp_coolant += f * power * dz / ( m_rate * c_p)
        
        # Initialize the temperature profile
        T_radial = [Temp_coolant] # Temperature of the coolant (ideally at r = infinity)

        r_coolant_cladding = geom_data.cladding_outer_diameter[i_h] / 2
        r_cladding_gap = geom_data.cladding_outer_diameter[i_h] / 2 - geom_data.thickness_cladding[i_h]
        r_gap_fuel = geom_data.fuel_outer_diameter[i_h] / 2
        r_fuel_in = geom_data.fuel_inner_diameter[i_h] / 2

        # Index of the interfaces
        idx_fuel = np.argmin(np.abs(r_values - r_gap_fuel)) 
        idx_gap = np.argmin(np.abs(r_values - r_cladding_gap)) 
        idx_cladding = np.argmin(np.abs(r_values - r_coolant_cladding))

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
                    fuel.Percent_of_Theoretical_Density = 98
                    
                # Equiaxed Region
                elif T_radial[j-1] > 1600:
                    fuel.Percent_of_Theoretical_Density = 96

                # Thermal resistance of the fuel
                th_res = thermal_resistance_fuel(Burnup, fuel, T_radial[j-1]) 
                # Compute the temperature
                T_value = T_radial[idx_fuel] + power * th_res * (1 - (r - r_fuel_in)**2/(r_gap_fuel - r_fuel_in)**2)
            # In the gap
            elif r < r_cladding_gap and idx_fuel != idx_gap:
                # Thermal resistance of the gap
                th_res = thermal_resistance_gap(geom_data, gap, fuel, cladding, T_radial[j-1], He_percentage, T_fuel_out, i_h)
                # Compute the temperature
                T_value = T_radial[idx_gap] + power * th_res * np.log(r_cladding_gap / r)
                
            # In the cladding
            elif r < r_coolant_cladding:
                # Thermal resistance of the cladding
                th_res = thermal_resistance_cladding(geom_data, cladding, T_radial[j-1], i_h)
                # Compute the temperature
                T_value = T_radial[idx_cladding] + power * th_res * (dr / (r_coolant_cladding - r_cladding_gap))
                
            # In the coolant
            else:
                # Thermal resistance of the coolant
                th_res, velocity = coolant_thermohydraulics(geom_data, thermo_hyd_spec, coolant, T_radial[0], i_h)
                # Compute the temperature
                T_value = T_radial[j-1] + power * th_res * (dr / (r_values[0] - r_coolant_cladding))

            # Compute new value of T
            T_radial.append(T_value)
        Z.append(T_radial)
        velocity_vector.append(velocity)
    Z = np.array(Z)
    T_map = Temperature_Map(X, Y, Z)
    return T_map, velocity_vector

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
    cladding_outer_diameter_swelled = []
    thickness_cladding_swelled = []

    for idx_h, h in enumerate(T_map.h[:, 0]):
        r_coolant_cladding = geom_data.cladding_outer_diameter[idx_h] / 2
        r_cladding_gap = r_coolant_cladding - geom_data.thickness_cladding[idx_h]
        r_gap_fuel = geom_data.fuel_outer_diameter[idx_h] / 2

        # Check gap size
        current_gap_size = r_cladding_gap - r_gap_fuel

        # Compute the average temperature in the cladding region
        idx_start = np.argmin(np.abs(T_map.r[idx_h, :] - r_coolant_cladding))
        idx_end = np.argmin(np.abs(T_map.r[idx_h, :] - r_cladding_gap))
        r_vals = T_map.r[idx_h, idx_start:idx_end]

        # Compute average temperature in this region
        if len(r_vals) == 0:
            temperature_avg = get_temperature_at_point(h, r_cladding_gap, T_map)
        else:
            temperature = [get_temperature_at_point(h, r, T_map) for r in r_vals]
            temperature_avg = np.mean(temperature)

        # Compute void swelling based on the temperature and neutron flux
        phi = power_profile(h, thermo_hyd_spec, value='neutron_flux') * thermo_hyd_spec.uptime
        void_swelling_fraction = (
            1.5e-3 * np.exp(-2.5 * ((temperature_avg - 273 - 450) / 100) ** 2) * (phi / 1e22) ** 2.75
        )
        if void_swelling_fraction < 0:
            void_swelling_fraction = 0  # Ensure swelling is non-negative
        Volume_expansion_fission_gas.append(void_swelling_fraction)

        # Compute new cladding geometry due to swelling
        outer_r_before = geom_data.cladding_outer_diameter[idx_h] / 2
        outer_r_swelled = outer_r_before * (1 + 1 / 3 * void_swelling_fraction / 100)
        
        inner_r_before = outer_r_before - geom_data.thickness_cladding[idx_h]

        available_gap = inner_r_before - geom_data.fuel_outer_diameter[idx_h] / 2
        if available_gap > 0:
            inner_r_swelled = inner_r_before * (1 + 1 / 3 * void_swelling_fraction / 100)
        else:
            inner_r_swelled  = inner_r_before

        # Update swelled geometry
        cladding_outer_diameter_swelled.append(2 * outer_r_swelled)
        thickness_cladding_swelled.append(outer_r_swelled - inner_r_swelled)

    return Volume_expansion_fission_gas, cladding_outer_diameter_swelled, thickness_cladding_swelled

##################################################
# Central Void Formation
##################################################
def get_R_void(Fuel_Proprieties, R_col, R_eq, Geometry_Data):
        
    R_void = []
    density_af = Fuel_Proprieties.Density
    density_columnar=Fuel_Proprieties.Theoretical_Density * (1 - Fuel_Proprieties.Porosity_Columnar)
    density_equiaxed= Fuel_Proprieties.Theoretical_Density * (1 - Fuel_Proprieties.Porosity_Equiaxed)
    
    for r in range(len(R_col)):
        r_fuel_out = Geometry_Data.fuel_outer_diameter[r] / 2
        if R_col[r] > r_fuel_out:
            R_col[r] = r_fuel_out
        if R_eq[r] > r_fuel_out:
            R_eq[r] = r_fuel_out

        R_voidd = np.sqrt(R_col[r]**2 - R_eq[r]**2*(density_af/density_columnar) +(R_eq[r]**2 - R_col[r]**2)*(density_equiaxed/density_af))
        R_void.append(R_voidd)
        
    return R_void

##################################################
# Thermal Expansion
##################################################
def thermal_expansion(fuel, cladding, cold_geometrical_data, T_map):
    """
    Computes the thermal expansion of fuel and cladding along the height of the system.

    Parameters:
    - fuel: Object containing fuel properties, including Thermal_Expansion_Coeff.
    - cladding: Object containing cladding properties, including Thermal_Expansion_Coeff.
    - cold_geometrical_data: Object containing cold geometry data such as diameters and cladding thickness.
    - T_map: A map or function to retrieve temperature at a given point.
    - geom_data: Additional geometry-related data.

    Returns:
    - fuel_outer_hot: List of expanded outer diameters for the fuel.
    - fuel_inner_hot: List of expanded inner diameters for the fuel.
    - cladding_outer_hot: List of expanded outer diameters for the cladding.
    - thickness_cladding: List of updated cladding thickness values.
    """
    # Constants
    T_0 = 298.15  # Reference temperature (25 °C in Kelvin)
    
    # Initialize results
    fuel_outer_hot = []
    fuel_inner_hot = []
    cladding_outer_hot = []
    thickness_cladding = []

    # Expansion coefficients
    fuel_alpha = fuel.Thermal_Expansion_Coeff
    cladding_alpha = cladding.Thermal_Expansion_Coeff

    # Iterate through height values
    for i_h, h in enumerate(cold_geometrical_data.h_values):
        # Starting Dimensions
        fuel_outer = cold_geometrical_data.fuel_outer_diameter[i_h] / 2
        fuel_inner = cold_geometrical_data.fuel_inner_diameter[i_h] / 2
        cladding_outer = cold_geometrical_data.cladding_outer_diameter[i_h] / 2
        cladding_inner = cladding_outer - cold_geometrical_data.thickness_cladding[i_h]
        avg_radius_cladding = (cladding_outer + cladding_inner) / 2

        # Temperature at average radii
        T_fuel_0 = get_temperature_at_point(h, fuel_inner, T_map)
        T_fuel_E = get_temperature_at_point(h, fuel_outer, T_map)
        T_fuel = (T_fuel_0 + T_fuel_E) / 2
        T_cladding = get_temperature_at_point(h, avg_radius_cladding, T_map)

        # Calculate expanded fuel inner radius
        fuel_inner_expanded = fuel_inner * (1 + fuel_alpha * (T_fuel - T_0))
        fuel_inner_hot.append(2 * fuel_inner_expanded)

        # Calculate expanded fuel dimensions
        fuel_outer_expanded = fuel_outer * (1 + fuel_alpha * (T_fuel - T_0))
        fuel_outer_hot.append(2 * fuel_outer_expanded)

        # Calculate expanded cladding outer radius
        cladding_outer_expanded = cladding_outer * (1 + cladding_alpha(T_cladding))
        cladding_outer_hot.append(2 * cladding_outer_expanded)

        # Calculate expanded cladding inner radius
        cladding_inner_expanded = cladding_inner * (1 + cladding_alpha(T_cladding))

        # Update possible conflicting dimensions:
        # Fuel outer diameter and Cladding inner diameter
        available_gap = cladding_inner_expanded - fuel_outer_expanded
        if available_gap <= 0:
            max_thickness = cladding_outer_expanded - fuel_outer_expanded
            thickness_cladding.append(max_thickness)
        else:
            # Compute cladding thickness normally
            thickness_cladding.append(cladding_outer_expanded - cladding_inner_expanded)
        
    return fuel_outer_hot, fuel_inner_hot, cladding_outer_hot, thickness_cladding

##################################################
# Fission Gas Production
##################################################
# Size of bubble radius
a = 10e-6  # m

# Fission yield
fission_yield = 0.3

# Diffusivity evaluation parameters [Matzke, 1980]
d_0 = 5e-8  # m^2/s
q = 40262

def fission_gas_production(h_plenum, Fuel_Proprieties, ThermoHydraulics, Geometrical_Data, T_map):
    d_0 = 5e-8  # m^2/s
    q = 40262
    """
    Function to compute the rate theory fission gas calculations and generate relevant plots.
    Outputs:
        He_percentage: Percentage of helium present inside the fuel still trapped in
        new_p_gas: Gas pressure in the plenum
    """
    # Calculate molar mass
    molar_mass = nf.mixture(Fuel_Proprieties.Molar_Mass, Fuel_Proprieties.Qualities, normalization_cond='normalize')  # g/mol

    # Calculate macroscopic cross sections
    macro_235 = nf.macro(Fuel_Proprieties.Micro_Fission[0], Fuel_Proprieties.Density, molar_mass)  # cm^-1
    macro_238 = nf.macro(Fuel_Proprieties.Micro_Fission[1], Fuel_Proprieties.Density, molar_mass)  # cm^-1
    macro_pu = nf.macro(Fuel_Proprieties.Micro_Fission[3], Fuel_Proprieties.Density, molar_mass)  # cm^-1

    fission_xs = nf.mixture([macro_235, macro_238, 0, macro_pu], Fuel_Proprieties.Qualities, normalization_cond='normalize')  # cm^-1

    # Compute average neutron flux
    # power profile
    q_values = [power_profile(h, ThermoHydraulics) for h in Geometrical_Data.h_values]
    peak_to_average_ratio = max(q_values) / np.mean(q_values)
    average_neutron_flux = ThermoHydraulics.neutron_flux_peak / peak_to_average_ratio

    # Calculate average fission rate
    avg_fission_rate = average_neutron_flux * fission_xs * 1e6  # [fissions/m^3 s]

    # Diffusivity function
    diffusivity = lambda temperature: d_0 * np.exp(-q / temperature)  # m^2/s

    # Compute the average maximum fuel temperature
    fuel_inner_diameter_avg = np.mean(Geometrical_Data.fuel_inner_diameter)
    temperature_max_average_fuel = np.mean(
        [get_temperature_at_point(point, fuel_inner_diameter_avg / 2, T_map)
        for point in Geometrical_Data.h_values]
    )

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
    total_fission_gas = P_lambda(ThermoHydraulics.uptime)  # Total amount of fission gas produced

    # Total amount of fission gas inside the grains
    n_grains_pellet = 1e5  # Number of grains in a pellet
    n_pellets_pin = round(850e-3 / Geometrical_Data.fuel_pellet_height)  # Number of pellets in a pin
    total_fission_gas_grains = GM_final * n_grains_pellet * n_pellets_pin

    # Total amount of fission gas released in the plenum
    total_fission_gas_released = total_fission_gas - total_fission_gas_grains

    # Calculate the percentage of helium trapped inside the fuel
    He_percentage = (total_fission_gas_grains / total_fission_gas)

    # Corresponding volume to accommodate gases
    r_cladding_gap = np.mean(np.array(Geometrical_Data.cladding_outer_diameter) / 2 - np.array(Geometrical_Data.thickness_cladding))
    r_gap_fuel = np.mean(np.array(Geometrical_Data.fuel_outer_diameter) / 2)
    
    V_plenum = (np.pi * r_cladding_gap**2 * h_plenum) + (np.pi * (r_cladding_gap**2 - r_gap_fuel**2) * Geometrical_Data.h_values[-1])

    # Find initial quantity of He present in the plenum
    initial_moles_he = Geometrical_Data.Initial_Gas_Pressure * V_plenum / (8.314 * (Geometrical_Data.Initial_Gas_Temperature))  # moles

    # Find the additional moles of fission gases released in the plenum
    fuel_outer_diameter_avg = np.mean(Geometrical_Data.fuel_outer_diameter)
    V_pin = np.pi * (fuel_outer_diameter_avg / 2)**2 * Geometrical_Data.h_values[-1]  # m^3
    additional_moles_fg = (total_fission_gas_released * V_pin) / 6.022e23  # moles

    # Find the total moles of gases in the plenum
    total_moles_gas = initial_moles_he + additional_moles_fg  # moles

    # Find the new pressure in the plenum
    new_p_gas = total_moles_gas * 8.314 * (Geometrical_Data.Initial_Gas_Temperature) / V_plenum  # Pa

    return He_percentage, new_p_gas

##################################################
# Main Function
##################################################
def initialize_params(delta):
    """
    Initialization function to set up parameters.
    
    Returns:
        dict: Dictionary containing initialized parameters.
    """
    # Material: Cladding
    # 15-15, Ti stabilized, cold worked stainless steel
    Cladding_Proprieties = Material_Proprieties(
        Elements=["Cr", "Ni", "Mo", "Mn", "Si", "Ti", "C", "B [ppm]"],
        Qualities=[15, 15, 1.5, 1.5, 0.9, 0.4, 0.09, 60],
        Density=lambda eth: 7900 * (1 + eth)**-3,  # kg/m^3
        Thermal_Conductivity=lambda t: 13.95 + 0.01163 * t,  # W/m K
        Emissivity = 0.32, # -
        Thermal_Expansion_Coeff=lambda t: -3.101e-4 + 1.545e-5 * (t-273.15) + 2.75e-9 * (t-273.15)**2,  # 1/°C --> 1/K
        Specific_Heat=500,  # Approximate value in J/kg K for steel
        Melting_Temperature=1673,  # K
        
        Youngs_Modulus=lambda t: 202.7 - 0.08167 * t,  # GPa
        Poissons_Ratio=lambda t: 0.277 + 6e-5 * t,  # Dimensionless
        Yield_Stress=lambda t: 555.5 - 0.25 * t if t < 600 else (405.5 - 0.775 * (t - 600) if t < 1000 else 345.5 - 0.25 * t),  # MPa
        Ultimate_Tensile_Strength=lambda t: 700 - 0.3125 * t if t < 600 else (512.5 - 0.969 * (t - 600) if t < 1000 else 437.5 - 0.3125 * t),  # MPa
    )

    # Material: Fuel
    # Homogeneous MOX fuel

    # Qualities have to be converted from weight to molar (won't change in hot condition):
    mol_qual = nf.w2mol([0.711, 0.29], [235 + 2*16, 239 + 2*16])  # UO2, PuO2

    Fuel_Proprieties = Material_Proprieties(
        Fission_Yield=0.3,
        Elements=["U-235", "U-238", "O-16", "Pu"],
        Qualities=[mol_qual[0] * (1 - mol_qual[1]), (1 - mol_qual[0]) * (1 - mol_qual[1]), 2, mol_qual[1]], # Molar fractions
        Micro_Fission = [1.047756375, 0.55801001, 0, 1.689844625],  # barn
        Theoretical_Density=11.31, # g/cm^3
        Percent_of_Theoretical_Density = 94.5, # %
        Porosity_Columnar = 0.02,
        Porosity_Equiaxed = 0.04,
        Molar_Mass=[235, 238, 16, 239],  # g/mol
        Thermal_Conductivity=lambda k_inf, beta: 1.755 + (k_inf - 1.755) * math.exp(-beta),  # W/m K
        Emissivity = 0.79, # -
        Thermal_Expansion_Coeff=1.2e-5,  # 1/°C
        Specific_Heat=270,  # Approximate value in J/kg K for MOX fuel
        Melting_Temperature=lambda pu, x, beta: 2964.92 + ((3147 - 364.85 * pu - 1014.15 * x) - 2964.92) * math.exp(-beta/41.01),  # K
        Oxigen_to_metal_ratio = 1.957, # -
        Grain_diameter = 10 * 1e-6,  # m
        Youngs_Modulus=lambda t, p: (22.43 * 10**4 - 31.19 * t) * (1 - 2.6 * p),  # MPa
        Poissons_Ratio=0.32,  # Dimensionless
    )

    # Material: Coolant (Sodium)
    Coolant_Proprieties = Material_Proprieties(
        Elements=["Na"],
        Qualities=[1],
        Density=lambda t: 954.1579 + ((t-273) * 9/5 +32) *( ((t-273) * 9/5 +32) * ((((t-273) * 9/5 +32) * 0.9667e-9 - 0.46e-5)) - 0.1273534),  # kg/m^3 (t * 9/5 +32) is the convertion from C to K to F
        Viscosity=lambda t: (np.exp(813.9 / t -2.530 ))/1000,  # Pa s
        Thermal_Conductivity=lambda t: 110 - 0.0648 * t + 1.16e-5 * t**2,  # W/m K
        Specific_Heat=lambda t: 1608 - 0.7481 * t + 3.929e-4 * t**2,  # J/kg K
        
        Melting_Temperature=98,  # °C
        Boiling_Temperature=882,  # °C
        Nusselt_Number=lambda pe: 7 + 0.025 * pe**0.8
    )

    # Material: Helium (Filling Gas)
    Helium_Proprieties = Material_Proprieties(
        Elements=["He"],
        Qualities=[1],
        Density=0.1786,  # kg/m^3 at STP
        Thermal_Conductivity=lambda t, x: (15.8e-4 * t**0.79)**x * (0.935e-4 * t**0.79)**(1-x),  # W/m K
        Specific_Heat=5193,  # J/kg K at constant pressure
        Thermal_Expansion_Coeff=3.66e-3  # Approximate value for helium in 1/°C
    )

    len_h = 10
    Geometrical_Data = GeometryData(
        Initial_Gas_Temperature=20 + 273,  # K
        Initial_Gas_Pressure=1e5,  # Pa
        fuel_outer_diameter=[5.42 * 1e-3]*len_h,  # m - GIVEN
        fuel_inner_diameter=[0.00 * 1e-3]*len_h,  # m
        cladding_outer_diameter=[6.55 * 1e-3]*len_h,  # m - GIVEN
        thickness_cladding=[delta]*len_h, # m
        pin_pitch=8.275 * 1e-3,  # m
        h_values = np.linspace(0, 0.85, len_h), # m
        fuel_pellet_height = 7e-3, # m
        fuel_roughness = 2e-6, # m
        cladding_roughness = 1e-6 # m
    )

    # Example of initializing Thermo-Hydraulic specifications
    heights_of_slice_centre = [42.5, 127.5, 212.5, 297.5, 382.5, 467.5, 552.5, 637.5, 722.5, 807.5] # mm
    ThermoHydraulics = ThermoHydraulicSpecs(
        coolant_inlet_temp=395 + 273,  # K
        coolant_inlet_pressure=1e5,  # Pa
        coolant_mass_flow_rate=0.049,  # kg/s
        q_linear_avg = 38.7e3,  #W/m,
        uptime = 360 * 24 * 3600 * 1,  # s
        h_peak_factor = [h * 1e-3 for h in heights_of_slice_centre],  # m
        peak_factors = [0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.498],
        neutron_flux_peak = 6.1e15  # Neutron Flux (> 100 keV) (10^15 n cm^-2 s^-1) at Peak Power Node    
    )

    params = {
        "Cladding_Proprieties": copy.deepcopy(Cladding_Proprieties),
        "Fuel_Proprieties": copy.deepcopy(Fuel_Proprieties),
        "Coolant_Proprieties": copy.deepcopy(Coolant_Proprieties),
        "Helium_Proprieties": copy.deepcopy(Helium_Proprieties),
        "Geometrical_Data_Cold": copy.deepcopy(Geometrical_Data),
        "Geometrical_Data": copy.deepcopy(Geometrical_Data),
        "ThermoHydraulics": copy.deepcopy(ThermoHydraulics),
    }
    #delete old vars
    del Cladding_Proprieties, Fuel_Proprieties, Coolant_Proprieties, Helium_Proprieties, Geometrical_Data, ThermoHydraulics
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
    HM_molar_mass = nf.mixture(HM_molar_mass, HM_qualities, normalization_cond='normalize') * 1e-3 # g/mol --> kg/mol
    HM_mass = MOX_mols * HM_molar_mass  # kg of Heavy Metal

    Burnup = ThermoHydraulics.uptime / (24 * 3600) \
        * ThermoHydraulics.q_linear_avg * Geometrical_Data_Cold.h_values[-1] \
        / HM_mass
    Burnup = Burnup * 1e-6  # Wd/kgU * 1e-9/1e-3 --> GWd/t(HM)
    return Burnup

def check_gap_closure(Geometrical_Data):
    Closure = False
    for i_h in range(len(Geometrical_Data.h_values)):
        r_cladding_gap = Geometrical_Data.cladding_outer_diameter[i_h] / 2 - Geometrical_Data.thickness_cladding[i_h]
        r_gap_fuel = Geometrical_Data.fuel_outer_diameter[i_h] / 2
        if r_gap_fuel >= r_cladding_gap:
            Closure = True
    return Closure

def mechanical_analysis(ThermoHydraulics, Cladding_Proprieties, Geometrical_Data, T_map, plenum_pressure):
    # Temperatures on the inner cladding 
    T_cladding_inner = [get_temperature_at_point(h, r, T_map) for h, r in zip(Geometrical_Data.h_values, Geometrical_Data.cladding_outer_diameter)]
    plastic_strain = []
    time_to_rupture = []
    for i, T_ref in enumerate(T_cladding_inner):
        Yield_stress = Cladding_Proprieties.Yield_Stress(T_ref-273.15)
        Ultimate_Tensile_Strength = Cladding_Proprieties.Ultimate_Tensile_Strength(T_ref-273.15)
        
        # Cladding radii
        r_cladding_inner = Geometrical_Data.cladding_outer_diameter[i] / 2 - Geometrical_Data.thickness_cladding[i]
        r_cladding_mid = Geometrical_Data.cladding_outer_diameter[i] / 2 - Geometrical_Data.thickness_cladding[i] / 2
        r_cladding_outer = Geometrical_Data.cladding_outer_diameter[i] / 2
        
        # Pressure inside the gap
        pressure_in = plenum_pressure
        # Pressure in the coolant
        pressure_out = ThermoHydraulics.coolant_inlet_pressure
        
        flag = True # Set flag for plastic strain on True by default
        
        #####################################
        # Mariotte Criterion
        Hoop_stress = r_cladding_mid * (pressure_in - pressure_out) * 1e-6 / Geometrical_Data.thickness_cladding[i] 
        if Hoop_stress < Yield_stress:
            flag = False
        #####################################
        
        #####################################
        # Lamè Criterion
        Costant_1 = 2 * (pressure_out - pressure_in) * r_cladding_inner**2 * r_cladding_outer**2 / (r_cladding_outer**2 - r_cladding_inner**2)
        Costant_2 = pressure_in + Costant_1 / (2 * r_cladding_inner**2)
        # Expressions obtained from stress analysis in Verification.ipynb
        stress_r = -Costant_1/(2*r_cladding_mid**2) + Costant_2
        sigma_theta = Costant_1/(2*r_cladding_mid**2) + Costant_2
        sigma_z = 2*Costant_2
        # Take the maximum difference between couples of the three stresses
        Lame_stress = max([abs(stress_r - sigma_theta), abs(stress_r - sigma_z), abs(sigma_theta - sigma_z)]) 
        Lame_stress = Lame_stress * 1e-6 # MPa
        if Lame_stress < 2/3*Yield_stress and Lame_stress < 1/3*Ultimate_Tensile_Strength:
            flag = False
        #####################################
        
        # Take as reference stress the maximum between the two criteria
        Stress = max(Hoop_stress, Lame_stress)
        
        #####################################
        # Evaluate effects:
        # Plastic Strain
        # Time to rupture due to creep
        #
        # Note: The two criteria are not exclusive, the plastic strain is computed if one of the two is satisfied
        # We have observed that the Mariotte stress is slightly higher than the Lamè (by decimals) 
        # but the limits on the Lamè criterion are much more restrictive (~33% more restrictive)
        if flag:
            plastic_strain.append(10) # Placeholder value to trigger the warning
            time_to_rupture.append(0) # No evaluation of the time to rupture since we are in the plastic regime
        else:
            # No plastic strain --> Evaluate rupture due to creep
            LM_parameter = (2060 - Stress)/0.095 # Larson-Miller Parameter with the equivalent stress

            temp = 10**((LM_parameter / T_ref) - 17.125) # hours
            temp = temp / (365*24) # convert to years
            time_to_rupture.append(temp)

            plastic_strain.append(0) # No plastic strain
        #####################################
        
    return plastic_strain, time_to_rupture

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
    R_void = get_R_void(params["Fuel_Proprieties"], R_columnar, R_equiaxied, Geometrical_Data)
    # Consider only void if it is enlarging
    for i_r, r in enumerate(R_void):
        if r > Geometrical_Data.fuel_inner_diameter[i_r] / 2:
            Geometrical_Data.fuel_inner_diameter = [2 * r for r in R_void]

    # Void Swelling
    void_swell, \
    Geometrical_Data.cladding_outer_diameter, \
    Geometrical_Data.thickness_cladding = void_swelling(previous_T_map, params["Geometrical_Data_Cold"], params["ThermoHydraulics"])

    # Fission Gas Production
    He_percentage, Gas_Pressure = fission_gas_production(h_plenum, params["Fuel_Proprieties"], params["ThermoHydraulics"], Geometrical_Data, previous_T_map)

    # Plastic Strain
    Plastic_Strain, Time_To_Rupture = mechanical_analysis(params["ThermoHydraulics"], params["Cladding_Proprieties"], Geometrical_Data, previous_T_map, Gas_Pressure)

    # Temperature Map
    T_map, Coolant_Velocity = temperature_map(params["Coolant_Proprieties"], params["Cladding_Proprieties"], params["Helium_Proprieties"], 
                              params["Fuel_Proprieties"], params["ThermoHydraulics"], Geometrical_Data, T_fuel_out, Burnup, He_percentage)
    idx_fuel = np.argmin(np.abs(T_map.r[5, :] - Geometrical_Data.fuel_outer_diameter[0]/2))
    T_fuel_out = T_map.T[5, idx_fuel]

    return T_map, T_fuel_out, Geometrical_Data, He_percentage, Gas_Pressure, Coolant_Velocity, void_swell, Plastic_Strain, Time_To_Rupture



    ############################################################################

def main(delta, h_plenum, settings):
    params = initialize_params(delta)
    #print("Parameters initialized.")

    Burnup = compute_burnup(params)
    #print("Burnup:", Burnup)

    # Initial values
    T_fuel_out = settings['T_fuel_out_initial']  # K (Initial guess)
    He_percentage = 1  # Initial Value

    # Set cladding thickness
    Geometrical_Data = params["Geometrical_Data"]

    residual = 1 # Placeholder for residual
    j = 0

    # Initial temperature map computed on cold geometrical data
    previous_T_map, _ = temperature_map(params["Coolant_Proprieties"], params["Cladding_Proprieties"], params["Helium_Proprieties"], 
                            params["Fuel_Proprieties"], params["ThermoHydraulics"], params["Geometrical_Data"], T_fuel_out, Burnup, He_percentage)
    ############################################################################
    # ITERATIVE LOOP
    if settings["hot_run"]:
        Closure = False
        while residual > settings['residual_threshold'] and  Closure == False and j<settings['run_limiter']:
            j += 1

            T_map, T_fuel_out, Geometrical_Data, He_percentage, Plenum_Pressure, Coolant_Velocity, Void_Swelling, Plastic_Strain, Time_to_Rupture = update_temperatures(params, Geometrical_Data, T_fuel_out, Burnup, He_percentage, h_plenum, previous_T_map)

            residual = np.mean(np.abs(T_map.T - previous_T_map.T)) / np.mean(previous_T_map.T)
            previous_T_map = copy.deepcopy(T_map)

            # gap_thickness = np.array(Geometrical_Data.cladding_outer_diameter) / 2 - np.array(Geometrical_Data.thickness_cladding) - np.array(Geometrical_Data.fuel_outer_diameter) / 2
            # print("Iteration:", j, "Residual:", residual, "Gap Thickness:", gap_thickness)

            Closure = check_gap_closure(Geometrical_Data)
            #if Closure:
                #print("Gap Closure Detected")

    else:
        T_map = previous_T_map
        print("Hot run disabled, disabling print results, enabling plotting")

    return T_map, Geometrical_Data, He_percentage, Plenum_Pressure, Coolant_Velocity, Void_Swelling, params, Burnup, coolant_infinity_limit, Burnup, Plastic_Strain, Time_to_Rupture