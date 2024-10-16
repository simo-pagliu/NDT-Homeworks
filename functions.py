#include libraries
import numpy as np

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
