#include libraries
import numpy as np
import nuclei_func as nf

# Define the class for the Material Properties
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
def power_profile_step(h, h_from_c, peak_factors, q_linear_avg):
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
    return power_values[-1]  # Return the last power value for h = 850

##################################################
# Radial Temperature Profile
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
    beta = 1  # You might need to provide the actual value of beta

    # Calculate k_0
    k_0 = (1 / (A + B * temperature) + D / (temperature**2) * np.exp(-E / temperature)) * (1 - porosity)**2.5

    # Calculate final thermal conductivity k
    k = 1.755 + (k_0 - 1.755) * np.exp(-Burnup / 128.75)

    # Calculate the thermal resistance of the fuel
    thermal_resistance = 1 / (4 * np.pi * k)
    return thermal_resistance