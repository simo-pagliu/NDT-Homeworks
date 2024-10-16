#include libraries
import numpy as np
import nuclei_func as nf

# Function to calculate the hydraulic flow parameters
def hydraulic_flow(mass_flow_rate, density, pitch, diameter_out):

    passage_area = ( 1/2 * pitch * pitch * np.sin(np.pi/3) ) - (1/2 * np.pi * (diameter_out/2)**2)

    # Calculate the velocity of the fluid
    velocity = mass_flow_rate / (density * passage_area)

    return velocity, passage_area


def hydraulic_diameter(passage_area, pitch, radius):
    
    # Calculate the hydraulic diameter
    wetted_perimeter = 3 * pitch - 6 * radius + np.pi * radius

    hydraulic_diameter = (4 * passage_area) / (wetted_perimeter)

    return hydraulic_diameter

def heat_trans_coefficient(density, diameter_out, mass_flow_rate, pitch, radius, viscosity, thermal_conductivity, c_p):

    velocity, passage_area= hydraulic_flow(mass_flow_rate, density, pitch, diameter_out)
    d_h = hydraulic_diameter(passage_area, pitch, radius)

    # Parameters
    reynolds = (density * velocity * d_h) / viscosity
    prandtl = c_p * viscosity / thermal_conductivity
    peclet = reynolds * prandtl
    nusselt = 7 + 0.025 *peclet**0.8

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

def fuel_mixture(micros, qualities, densities, molar_masses):

    macros = []
    for i in range(len(micros)):
        macros.append(nf.macro(micros[i], densities[i], molar_masses[i]))

    # Calculate the mixture macroscopic cross section
    mixture_xs = nf.mixture(macros, qualities)

    return mixture_xs
        

    #fission cross section 

def peak_power(peak_flux, micros, qualities, densities, molar_masses, energy_per_fission, radius, active_length):

    volume = np.pi * radius**2 * active_length
    fission_xs = fuel_mixture(micros, qualities, densities, molar_masses)
    # Calculate the peak power
    peak_power = peak_flux * fission_xs * energy_per_fission * volume 
    peak_power_linear = peak_power / active_length

    return peak_power_linear