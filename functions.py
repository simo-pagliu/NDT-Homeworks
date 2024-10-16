#include libraries

import numpy as np

# Function to calculate the hydraulic flow parameters
def hydraulic_flow(mass_flow_rate, density, pitch, diameter_out):

    passage_area = ( 1/2 * pitch * pitch * np.sin(np.pi/3) ) - (1/2 * np.pi * (diameter_out/2)**2)

    # Calculate the velocity of the fluid
    velocity = mass_flow_rate / (density * passage_area)

    return velocity, passage_area
