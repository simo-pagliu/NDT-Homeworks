#include libraries
import numpy as np

# Define the class for the Material Properties
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
        self.nominal_size = (cladding_outer_diameter - fuel_outer_diameter)/2 - thickness_cladding
        self.effective_gap_size = self.nominal_size + fuel_roughness + cladding_roughness


# Thermo Hydraulics Specs
class ThermoHydraulicSpecs:
    def __init__(self, coolant_inlet_temp, coolant_inlet_pressure, coolant_mass_flow_rate, q_linear_avg, uptime, h_peak_factor, peak_factors, neutron_flux_peak):
        self.coolant_inlet_temp = coolant_inlet_temp  # K
        self.coolant_inlet_pressure = coolant_inlet_pressure  # Pa
        self.coolant_mass_flow_rate = coolant_mass_flow_rate  # kg/s
        self.q_linear_avg = q_linear_avg  # W/m
        self.uptime = uptime  # s
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
# Hydraulic Flow
##################################################
# Function to calculate the hydraulic flow parameters
def hydraulic_flow(thermo_hyd_spec, geom_data, coolant, temperature):
    passage_area = (1/2 * (geom_data.pin_pitch ** 2) * np.sin(np.pi/3)) - (1/2 * np.pi * (geom_data.cladding_outer_diameter/2)**2)

    # Calculate the velocity of the fluid
    velocity = thermo_hyd_spec.coolant_mass_flow_rate / (coolant.Density(temperature) * passage_area)

    # Calculate the hydraulic diameter
    wetted_perimeter = np.pi * geom_data.cladding_outer_diameter / 2

    hydraulic_diameter = (4 * passage_area) / (wetted_perimeter)

    return velocity, passage_area, hydraulic_diameter

def heat_trans_coefficient(geom_data, thermo_hyd_spec, coolant, temperature):
    # Import material properties
    density = coolant.Density(temperature)
    viscosity = coolant.Viscosity(temperature)
    thermal_conductivity = coolant.Thermal_Conductivity(temperature)
    c_p = coolant.Specific_Heat(temperature)

    # Calculate the velocity and passage area
    velocity, passage_area, d_h = hydraulic_flow(thermo_hyd_spec, geom_data, coolant, temperature)

    # Adimensional numbers
    reynolds = (density * velocity * d_h) / viscosity
    prandtl = c_p * viscosity / thermal_conductivity
    peclet = reynolds * prandtl
    nusselt = coolant.Nusselt_Number(peclet)

    # HTC calculation
    htc = nusselt * thermal_conductivity / d_h

    return htc, reynolds, prandtl, peclet, nusselt, d_h

# Define a function that outputs a power value for a given h
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
def thermal_resistance_coolant(geom_data, coolant_htc):
    radius_cladding_out = geom_data.cladding_outer_diameter / 2
    thermal_resistance = 1 / (2 * np.pi * radius_cladding_out * coolant_htc)
    return thermal_resistance

def thermal_resistance_cladding(geom_data, cladding, temperature):
    radius_cladding_out = geom_data.cladding_outer_diameter / 2
    radius_cladding_in = radius_cladding_out - geom_data.thickness_cladding
    k = cladding.Thermal_Conductivity(temperature)
    thermal_resistance = np.log(radius_cladding_out / radius_cladding_in) / (2 * np.pi * k)
    return thermal_resistance

def thermal_resistance_gap(geom_data, helium, fuel, cladding):
    radius_gap_in = geom_data.fuel_outer_diameter / 2

    # Conductive heat transfer
    conduction = lambda t_gas: helium.Thermal_Conductivity(t_gas) / geom_data.effective_gap_size

    # Radiative heat transfer
    radiation = lambda t_fuel_out: 4 * 5.67e-8 * (t_fuel_out**3) / (1/(fuel.Emissivity) + 1/(cladding.Emissivity) - 1)

    # Overall heat transfer coefficient
    htc = lambda t_gas, t_fuel_out: conduction(t_gas) + radiation(t_fuel_out)

    # Calculate the thermal resistance of the gap
    thermal_resistance = lambda t_gas, t_fuel_out: 1 / (2 * np.pi * radius_gap_in * htc(t_gas, t_fuel_out))

    return thermal_resistance

def thermal_resistance_fuel(Burnup, fuel):
    A = 0.01926 + 1.06e-6 * fuel.Oxigen_to_metal_ratio + 2.63e-8 * fuel.Molar_Mass[-1]
    B = 2.39e-4 + 1.37e-13 * fuel.Molar_Mass[-1]
    D = 5.27e9
    E = 17109.5

    # Calculate k_0
    k_0 = lambda temperature: (1 / (A + B * temperature) + D / (temperature**2) * np.exp(-E / temperature)) * (1 - fuel.Porosity)**2.5

    # Calculate final thermal conductivity k
    k = lambda temperature: 1.755 + (k_0(temperature) - 1.755) * np.exp(-Burnup / 128.75)

    # Calculate the thermal resistance of the fuel
    thermal_resistance = lambda temperature: 1 / (4 * np.pi * k(temperature))
    return thermal_resistance

##################################################
# Temperature Profiles
##################################################
def axial_T_profile_coolant(Temp_old, thermo_hyd_spec, power, h, dz, coolant):
    f = 1/2
    c_p = coolant.Specific_Heat(Temp_old)
    m_rate = thermo_hyd_spec.coolant_mass_flow_rate

    power = power_profile(h, thermo_hyd_spec)

    T_now = Temp_old + power * dz / ( m_rate / f * c_p)
    return T_now

def radial_temperature_profile(Temp_0, power, r_plot, geom_data, Resistances, T_fuel_out, Burnup):
    # Initialize the temperature profile
    T_radial = [Temp_0] # Temperature of the coolant (ideally at r = infinity)

    # Create Limits
    r_coolant_cladding = geom_data.cladding_outer_diameter / 2
    r_cladding_gap = geom_data.cladding_outer_diameter / 2 - geom_data.thickness_cladding
    r_gap_fuel = geom_data.fuel_outer_diameter / 2
    r_fuel_in = geom_data.fuel_inner_diameter / 2

    # Index of the interfaces
    idx_fuel = np.argmin(np.abs(r_plot - r_gap_fuel))
    idx_gap = np.argmin(np.abs(r_plot - r_cladding_gap))

    # Compute the temperature profile
    for j, r in enumerate(r_plot[1:], start=1): 
        dr = r_plot[j-1] - r_plot[j]

        # In the pellet
        if r < r_gap_fuel:
            th_res = Resistances.Fuel(T_radial[j-1])
            # Compute the temperature
            if r_fuel_in == 0:
                T_value = T_radial[idx_fuel] + power * th_res * (1 - (r/r_gap_fuel)**2)
            else:
                T_value = T_radial[idx_fuel] + power * th_res * np.log(r_gap_fuel / r)
        
        # In the gap
        elif r < r_cladding_gap:
            th_res = Resistances.Gap(T_radial[j-1], T_fuel_out)
            # Compute the temperature
            T_value = T_radial[idx_gap] + power * th_res * np.log(r_cladding_gap / r)
            
        # In the cladding
        elif r < r_coolant_cladding:
            th_res = Resistances.Cladding
            T_value = T_radial[j-1] + power * th_res * (dr / (r_coolant_cladding - r_cladding_gap))
            
        # In the coolant
        else:
            th_res = Resistances.Coolant
            T_value = T_radial[j-1] + power * th_res * (dr / (r_plot[0] - r_coolant_cladding))

        # Compute new value of T
        T_radial.append(T_value)
        
    return T_radial

def temperature_profile_3D(r_values, Resistances, coolant, thermo_hyd_spec, geom_data, T_fuel_out, Burnup):
    h_values = geom_data.h_values

    # Compute the height of each slice
    dz = h_values[1] - h_values[0]
    

    # Create the meshgrid
    X, Y = np.meshgrid(r_values, h_values)
    Z = []
    
    # Initialize the temperature
    Temp_coolant = thermo_hyd_spec.coolant_inlet_temp

    # Compute the temperature profile for each height step
    for h in h_values: 
        q = power_profile(h, thermo_hyd_spec)
        # Compute temperature profile
        Temp_coolant = axial_T_profile_coolant(Temp_coolant, thermo_hyd_spec, q, h, dz, coolant)
        T_plot = radial_temperature_profile(Temp_coolant, q, r_values, geom_data, Resistances, T_fuel_out, Burnup)
        Z.append(T_plot)

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
# Void Formation
##################################################
def void_swelling(T_map, geom_data, thermo_hyd_spec):
    Volume_expansion_fission_gas = []
    r_coolant_cladding = geom_data.cladding_outer_diameter / 2
    r_cladding_gap = geom_data.cladding_outer_diameter / 2 - geom_data.thickness_cladding
    idx_start = np.argmin(np.abs(T_map.r[0, :] - r_coolant_cladding))
    idx_end = np.argmin(np.abs(T_map.r[0, :] - r_cladding_gap))
    r_vals = T_map.r[0, idx_start:idx_end]
    
    for h_idx, h in enumerate(T_map.h[:, 0]):
        phi = power_profile(h, thermo_hyd_spec, value = 'neutron_flux') * thermo_hyd_spec.uptime

        temperature = [get_temperature_at_point(h, r, T_map) for r in r_vals]
        temperature_avg = np.mean(temperature)

        temp = 1.5e-3 * np.exp(-2.5 * ((temperature_avg - 273 - 450) / 100) ** 2) * (phi / 1e22) ** 2.75
        if temp < 0:
            temp = 0
        Volume_expansion_fission_gas.append(temp)
    return Volume_expansion_fission_gas

def get_R_void(Fuel_Properties, R_col, R_eq):
        
    R_void = []
    density_af = Fuel_Properties.Density
    density_columnar=Fuel_Properties.Theoretical_Density*0.98
    density_equiaxed= Fuel_Properties.Theoretical_Density*0.95
    
    for r in range(len(R_col)):
        
        R_voidd = np.sqrt(R_col[r]**2 - R_eq[r]**2*(density_af/density_columnar) +(R_eq[r]**2 - R_col[r]**2)*(density_equiaxed/density_af))
        R_void.append(R_voidd)
        
    return R_void
##################################################
# Thermal Expansion
##################################################
def cold_to_hot_fuel(Fuel, Geometrical_Data, vars, h_vals):

    alpha = Fuel.Thermal_Expansion_Coeff
    T_0 = 25 + 273.15

    R_hot = []
    R_init = []
    T_hot= []

    for h in h_vals:

        T_hott = get_temperature_at_point(h, Geometrical_Data.fuel_outer_diameter/2, vars.T_map)
        T_hot.append(T_hott)
        R_initt = Geometrical_Data.fuel_outer_diameter/2
        R_init.append(R_initt)
        R_hott =R_initt * (1 + alpha * (T_hott - T_0))
        R_hot.append(R_hott)
    
    return R_hot, R_init, T_hot

def cold_to_hot_clad (Cladding, Geometrical_Data, vars, h_vals, flag):
    
    if flag=='inner':
        diam = Geometrical_Data.cladding_outer_diameter/2 - Geometrical_Data.thickness_cladding
    elif flag=='outer':
        diam = Geometrical_Data.cladding_outer_diameter/2
        
    strain = []
    R_hot = []
    R_init = []
    T_hot= []

    for h in h_vals:

        T_hott = get_temperature_at_point(h, diam, vars.T_map) 
        T_hot.append(T_hott)
        
        strainn = Cladding.Thermal_Expansion_Coeff(T_hott)
        strain.append(strainn)
        
        R_initt = diam
        R_init.append(R_initt)
        
        R_hott =R_initt * (1 + strainn)
        R_hot.append(R_hott)
    
    return R_hot, R_init, T_hot

##################################################
# Fission Gas Release
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

import nuclei_func as nf
def fission_gas_production(h_plenum, Fuel_Proprieties, ThermoHydraulics, Geometrical_Data, T_map):
    """
    Function to compute the rate theory fission gas calculations and generate relevant plots.
    Outputs:
        He_percentage: Percentage of helium present inside the fuel still trapped in
        new_p_gas: Gas pressure in the plenum
    """
    # Calculate molar mass
    molar_mass = nf.mixture(Fuel_Proprieties.Molar_Mass, Fuel_Proprieties.Qualities)  # g/mol

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
    He_percentage = (total_fission_gas_grains / total_fission_gas) * 100

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

## THERMAL CREEP

def rupture_strain(T_map, cladding_inner_diam):
    
    rupture_strain = []
    
    for h in T_map.h[:,0]:
        
        temperature=get_temperature_at_point(h, cladding_inner_diam, T_map)
        
        if temperature < 500:
            rupture_strain.append(0.1e-5)
        else:
            rupture_strain.append(8 +4.74e-3*(temperature-500) + 6.2e-5*(temperature-500)**2)
    
    return rupture_strain, temperature

def LMP(T_map, useful_life, mat_const,rupture_strain,cladding_inner_diam):
    
    LMP = []
    
    for i_h, h in enumerate(T_map.h[:,0]):
        
        temperature=get_temperature_at_point(h, cladding_inner_diam, T_map)
        
        LMP.append(temperature*(mat_const-np.log10(rupture_strain[i_h]/useful_life)))
    
    return LMP, temperature

def get_stress_from_graph(max_temperature, useful_life):

    rup_strain = (8 + 4.74e-3*(max_temperature-500) + 6.2e-5*(max_temperature-500)**2)*1e-2
    max_LMP = max_temperature*(17.125-np.log10(rup_strain/useful_life))
    max_stress = 2060 - 0.095*max_LMP
    
    return max_stress, max_LMP, rup_strain