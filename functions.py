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
        

def peak_power(peak_flux, fuel, energy_per_fission, radius, active_length):
    # Calculate the volume of the fuel   
    volume = np.pi * radius**2 * active_length
    # Calculate the fission cross section
    fission_xs = fuel_mixture(fuel)
    # Calculate the peak power
    peak_power = peak_flux * fission_xs * energy_per_fission * volume 

    #obtained value is in [W]
    return peak_power, fission_xs

