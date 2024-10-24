#include libraries
import numpy as np
import nuclei_func as nf
import matplotlib.pyplot as plt

# Define a class for the geometry data of the reactor components
class GeometryData:
    def __init__(self, fuel_outer_diameter, fuel_inner_diameter, cladding_outer_diameter, thickness_cladding):
        self.fuel_outer_diameter = fuel_outer_diameter  # mm
        self.fuel_inner_diameter = fuel_inner_diameter  # mm (0 if solid fuel pellet)
        self.cladding_outer_diameter = cladding_outer_diameter  # mm
        self.thickness_cladding = thickness_cladding  # mm

    def get_fuel_radius(self):
        return self.fuel_outer_diameter / 2

    def get_cladding_radius(self):
        return self.cladding_outer_diameter / 2

    def get_effective_gap_size(self, fuel_roughness, cladding_roughness):
        nominal_size = self.cladding_outer_diameter - self.fuel_outer_diameter - self.thickness_cladding
        return nominal_size + fuel_roughness + cladding_roughness

# Define a class for Thermo-Hydraulic specifications
class ThermoHydraulicSpecs:
    def __init__(self, coolant_inlet_temp, coolant_inlet_pressure, coolant_mass_flow_rate, pin_pitch):
        self.coolant_inlet_temp = coolant_inlet_temp  # K
        self.coolant_inlet_pressure = coolant_inlet_pressure  # MPa
        self.coolant_mass_flow_rate = coolant_mass_flow_rate  # kg/s
        self.pin_pitch = pin_pitch  # mm

# Modify the Material_Proprieties class to include specific parameters related to coolant
class Material_Proprieties:
    def __init__(self, Elements='', Qualities='', Density='', Molar_Mass='',Micro_Fission = '', Micro_Absorption = '', Viscosity='', Thermal_Conductivity='', Specific_Heat='', Thermal_Expansion_Coeff='', Melting_Temperature='', Boiling_Temperature='', Youngs_Modulus='', Poissons_Ratio='', Yield_Stress='', Ultimate_Tensile_Strength='', Nusselt_Number=''):
        self.Elements = Elements
        self.Qualities = Qualities
        self.Density = Density
        self.Molar_Mass = Molar_Mass
        self.Micro_Fission = Micro_Fission
        self.Micro_Absorption = Micro_Absorption
        self.Viscosity = Viscosity
        self.Thermal_Conductivity = Thermal_Conductivity
        self.Specific_Heat = Specific_Heat
        self.Thermal_Expansion_Coeff = Thermal_Expansion_Coeff
        self.Melting_Temperature = Melting_Temperature
        self.Boiling_Temperature = Boiling_Temperature
        self.Youngs_Modulus = Youngs_Modulus
        self.Poissons_Ratio = Poissons_Ratio
        self.Yield_Stress = Yield_Stress
        self.Ultimate_Tensile_Strength = Ultimate_Tensile_Strength
        self.Nusselt_Number = Nusselt_Number

# Function to calculate the hydraulic flow parameters
def hydraulic_flow(mass_flow_rate, pitch, diameter_out, coolant, temperature):
    passage_area = ( 1/2 * pitch * pitch * np.sin(np.pi/3) ) - (1/2 * np.pi * (diameter_out/2)**2)

    # Calculate the velocity of the fluid
    velocity = mass_flow_rate / (coolant.Density(temperature) * passage_area)

    
    # Calculate the hydraulic diameter
    wetted_perimeter = np.pi* diameter_out/2

    hydraulic_diameter = (4 * passage_area) / (wetted_perimeter)

    return velocity, passage_area, hydraulic_diameter

def heat_trans_coefficient(diameter_out, mass_flow_rate, pitch, coolant, temperature):
    # Import material properties
    density = coolant.Density(temperature)
    viscosity = coolant.Viscosity(temperature)
    thermal_conductivity = coolant.Thermal_Conductivity(temperature)
    c_p = coolant.Specific_Heat(temperature)

    # Calculate the velocity and passage area
    velocity, passage_area, d_h = hydraulic_flow(mass_flow_rate, pitch, diameter_out, coolant, temperature)

    # Adimensional numbers
    reynolds = (density * velocity * d_h) / viscosity
    prandtl = c_p * viscosity / thermal_conductivity
    peclet = reynolds * prandtl
    nusselt = coolant.Nusselt_Number(peclet)


    #HTC calculation
    htc = nusselt * thermal_conductivity / d_h

    return htc, reynolds, prandtl, peclet, nusselt, d_h   

# Define a function that outputs a power value for a given h
def power_profile(h, h_from_c, peak_factors, q_linear_avg):
    # Computes peak power such that the average power is q_linear_avg
    peak_value = q_linear_avg * len(peak_factors) / sum(peak_factors)

    # Compute power values for each interval
    power_values = [peak_value * factor for factor in peak_factors]

    # Find the interval that h belongs to
    length_h_from_c = len(h_from_c)
    interval_boundaries = [0] + [(h_from_c[i] + h_from_c[i + 1]) / 2 for i in range(length_h_from_c - 1)] + [850]

    for i in range(length_h_from_c):
        if interval_boundaries[i] <= h <= interval_boundaries[i + 1]:
            return power_values[i]

##################################################
# Thermal Resistances
##################################################
def thermal_resistance_coolant(radius_cladding_out, coolant_htc):
    # Calculate the thermal resistance of the coolant
    thermal_resistance = 1 / (2 * np.pi * radius_cladding_out * coolant_htc)
    return thermal_resistance

def thermal_resistance_cladding(radius_cladding_out, radius_cladding_in, cladding, temperature):
    # Calculate the thermal resistance of the cladding
    k = cladding.Thermal_Conductivity(temperature)
    thermal_resistance = np.log(radius_cladding_out / radius_cladding_in) / (2 * np.pi * k)
    return thermal_resistance

def thermal_resistance_gap(effective_size, nominal_size, radius_gap_in, 
                           helium, fuel, cladding, 
                           temperature_gas, temperature_fuel_out, temperature_cladding_in,
                           contact_pressure, cladding_hardness, C,
                           emissivity_fuel, emissivity_cladding):

    # Three components of the thermal resistance
    # Conductive heat transfer
    conduction = helium.Thermal_Conductivity(temperature_gas) / effective_size 

    # Radiative heat transfer
    radiation = 4 * 5.67e-8 * (temperature_fuel_out**3) / (1/(emissivity_fuel) + 1/(emissivity_cladding) - 1)

    # Contact heat transfer
    # contact = C * (fuel.Thermal_Conductivity(temperature_fuel_out) * 
    #                cladding.Thermal_Conductivity(temperature_cladding_in) * 
    #                contact_pressure) / ((fuel.Thermal_Conductivity(temperature_fuel_out) + 
    #                cladding.Thermal_Conductivity(temperature_cladding_in))*
    #                (np.sqrt(nominal_size)*cladding_hardness))
    
    # Overall heat transfer coefficient
    htc = conduction + radiation #+ contact

    # Calculate the thermal resistance of the gap
    thermal_resistance = 1/(2 * np.pi * radius_gap_in * htc)
    return thermal_resistance

def thermal_resistance_fuel(Burnup, temperature, Oxigen_to_metal_ratio, Pu_concentration, porosity):
    # Define constants
    A = 0.01926 + 1.06e-6 * Oxigen_to_metal_ratio + 2.63e-8 * Pu_concentration
    B = 2.39e-4 + 1.37e-13 * Pu_concentration
    D = 5.27e9
    E = 17109.5

    # Calculate k_0
    k_0 = (1 / (A + B * temperature) + D / (temperature**2) * np.exp(-E / temperature)) * (1 - porosity)**2.5

    # Calculate final thermal conductivity k
    k = 1.755 + (k_0 - 1.755) * np.exp(-Burnup / 128.75)

    # Calculate the thermal resistance of the fuel
    thermal_resistance = 1 / (4 * np.pi * k)
    return thermal_resistance

##################################################
# Temperature Profiles
##################################################
def axial_T_profile_coolant(temperature_inlet, power, h, dz, coolant, mass_flow_rate, heights_of_slice_centre, peak_factors, q_linear_avg):
    # Computes axial temperature profile
    f = 1/2
    c_p = coolant.Specific_Heat(temperature_inlet)
    
    power = power_profile(h, heights_of_slice_centre, peak_factors, q_linear_avg)

    T_now = temperature_inlet + power * (dz) / (mass_flow_rate / f * c_p)
    
    return T_now
  
def temperature_profile(power, thermal_resistance, T_init):
    # Calculate final temperature given power and thermal resistance
    return T_init + power * thermal_resistance

def get_resistance(r, fuel_pellet_inner_diameter, fuel_pellet_outer_diameter, thickness_cladding, cladding_outer_diameter, Resistances):
    # Return appropriate thermal resistance based on the radius
    if r < fuel_pellet_outer_diameter:
        if fuel_pellet_inner_diameter == 0:
            res = Resistances.Fuel * (- r**2)
        else:
            res = Resistances.Fuel * np.log(fuel_pellet_outer_diameter*1e-3/2 / r) / np.log(fuel_pellet_inner_diameter / fuel_pellet_outer_diameter)
        return res
    elif r < fuel_pellet_outer_diameter + thickness_cladding:
        return Resistances.Gap
    elif r < cladding_outer_diameter:
        return Resistances.Cladding
    else:
        return Resistances.Coolant

def compute_temperature_profile(Temp_0, reference_power_density, r_plot, Resistances, fuel_pellet_inner_diameter, fuel_pellet_outer_diameter, thickness_cladding, cladding_outer_diameter):
    # Compute temperature profile for a given power density and radius plot
    T_plot = [Temp_0]
    for j, r in enumerate(r_plot[1:], start=1):
        R = get_resistance(r, fuel_pellet_inner_diameter ,fuel_pellet_outer_diameter, thickness_cladding, cladding_outer_diameter, Resistances)
        T_plot.append(temperature_profile(reference_power_density, R, T_plot[j-1]))
    return T_plot

def plot_3d_temperature_profile(h_values, r_plot, heights_of_slice_centre, peak_factors, q_linear_avg, Resistances, coolant, temperature_inlet, mass_flow_rate, fuel_pellet_inner_diameter, fuel_pellet_outer_diameter, thickness_cladding, cladding_outer_diameter):
    # Create a 3D plot for temperature profiles at different heights and store the temperature matrix
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    dz = h_values[1] - h_values[0]

    # Arrays for 3D plotting
    X, Y = np.meshgrid(r_plot, h_values)
    Z = []

    Temp_0 = temperature_inlet
    for h in h_values:       
        # Compute power profile at this height
        q = power_profile(h, heights_of_slice_centre, peak_factors, q_linear_avg)

        # Get the initial temperature based on the height
        Temp_0 = axial_T_profile_coolant(Temp_0, q, h, dz, coolant, mass_flow_rate, heights_of_slice_centre, peak_factors, q_linear_avg)
        
        # Compute temperature profile for this height
        T_plot = compute_temperature_profile(Temp_0, q, r_plot, Resistances, fuel_pellet_inner_diameter, fuel_pellet_outer_diameter, thickness_cladding, cladding_outer_diameter)
        Z.append(T_plot)

    Z = np.array(Z)

    # Plotting the surface
    ax.plot_surface(X, Y, Z, cmap='viridis')
    ax.set_xlabel('Radius (mm)')
    ax.set_ylabel('Height (mm)')
    ax.set_zlabel('Temperature (K)')
    ax.set_title("3D Temperature Profile vs. Radius and Height (q_values)")
    plt.show()

    return X, Y, Z

def get_temperature_at_point(h, r, X, Y, Z):
    # Get temperature at a specific height and radius using precomputed matrix Z
    # Find the closest indices in X and Y to the input r and h
    h_idx = np.argmin(np.abs(Y[:, 0] - h))
    r_idx = np.argmin(np.abs(X[0, :] - r))

    # Extract the corresponding temperature from Z
    return Z[h_idx, r_idx]