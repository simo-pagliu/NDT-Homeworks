#include libraries
import numpy as np
import nuclei_func as nf
import matplotlib.pyplot as plt

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

    return velocity, passage_area


def hydraulic_diameter(passage_area, pitch, radius):
    
    # Calculate the hydraulic diameter
    wetted_perimeter = 3 * pitch - 6 * radius + np.pi * radius

    hydraulic_diameter = (4 * passage_area) / (wetted_perimeter)

    return hydraulic_diameter

def heat_trans_coefficient(diameter_out, mass_flow_rate, pitch, coolant, temperature):
    # Import material properties
    density = coolant.Density(temperature)
    viscosity = coolant.Viscosity(temperature)
    thermal_conductivity = coolant.Thermal_Conductivity(temperature)
    c_p = coolant.Specific_Heat(temperature)

    # Calculate the velocity and passage area
    velocity, passage_area= hydraulic_flow(mass_flow_rate, pitch, diameter_out, coolant, temperature)
    d_h = hydraulic_diameter(passage_area, pitch, diameter_out/2)

    # Adimensional numbers
    reynolds = (density * velocity * d_h) / viscosity
    prandtl = c_p * viscosity / thermal_conductivity
    peclet = reynolds * prandtl
    nusselt = coolant.Nusselt_Number(peclet)


    #HTC calculation
    htc = nusselt * thermal_conductivity / d_h

    return htc

def power_profile(peak, nodes_center, amplitude, z_extrapolated):
    components = [
        lambda z, i=i: peak * amplitude[i] * np.cos(np.pi * (z - nodes_center[i]) / z_extrapolated)
        for i in range(len(nodes_center))
    ]

    # Sum the lambda functions
    power_profile = lambda z: sum(component(z) for component in components)
    return power_profile

def fuel_mixture(fuel):
    micros = fuel.Micro_Fission
    qualities = fuel.Qualities
    density = fuel.Density # The density of the mixture is given, kg/m3 = g/cm3
    molar_masses = fuel.Molar_Mass
    
    macros = []
    for i in range(len(micros)):
        macros.append(nf.macro(micros[i], density, molar_masses[i]))

    # Calculate the mixture macroscopic cross section
    mixture_xs = nf.mixture(macros, qualities)
    
    return mixture_xs
        

def axial_T_profile(temperature_inlet, coolant, mass_flow_rate):
    # Computes axial temperature profile

    f = 1/2
    c_p = coolant.Specific_Heat(temperature_inlet)
    z_in = 0

    power = [0, 27.7120681, 35.7059339, 42.05257887, 46.4128693, 48.44767151, 47.62406109, 44.18427641, 38.85503255, 31.87856785, 24.12694041]
    z = [0, 0.0425, 0.1275, 0.2125, 0.2975, 0.3825, 0.4675, 0.5525, 0.6375, 0.7225, 0.8075]
    
    T_prof = []

    for i in range(len(z)):
        if i == 0:
            T_now = temperature_inlet
            T_prof.append(T_now)
        else: 
            dz = z[i] - z [i-1]
            T_old = T_now
            T_now = T_old + power[i] *1e3* (dz) / (mass_flow_rate/f*c_p)
            T_prof.append(T_now)
    
    return z, T_prof
