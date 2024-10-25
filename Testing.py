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