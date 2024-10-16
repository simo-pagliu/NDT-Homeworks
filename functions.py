#include libraries
import numpy as np

# Function to calculate the hydraulic flow parameters
def hydraulic_flow(mass_flow_rate, density, pitch, diameter_out):

    passage_area = ( 1/2 * pitch * pitch * np.sin(np.pi/3) ) - (1/2 * np.pi * (diameter_out/2)**2)

    # Calculate the velocity of the fluid
    velocity = mass_flow_rate / (density * passage_area)

    return velocity, passage_area

def reynolds_number(velocity, diameter, density, viscosity):
    # Calculate the Reynolds number
    reynolds = (density * velocity * diameter) / viscosity

    return reynolds

def hydraulic_diameter(passage_area, pitch, radius):
    
    # Calculate the hydraulic diameter
    wetted_perimeter = 3 * pitch - 6 * radius + np.pi * radius

    hydraulic_diameter = (4 * passage_area) / (wetted_perimeter)

    return hydraulic_diameter