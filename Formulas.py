# Heat

# Heat Generation
q_i = m_flow * c_p * delta_T # [W/m]
q_i = m_flow * delta_h # [W/m]
q_ii = q_i / (2* math.pi * r) # Heat Generation [W/m2]

# Heat Transfer
q_ii = htc * delta_T # Heat transfer [W/m2]

# Thermal Expansion same formula valid for linear, radial (azimuthal really...) and volumetric expansion
delta_L = L_0 * alpha * delta_T # Linear Expansion [m]

# Pressure Drops
delta_p_f = f_darcy * (length / D_hydraulic) * (density * velocity**2) / 2 # Frictional Pressure Drop Pa
delta_p_g = (density * g * length) # Gravitational Pressure Drop Pa
delta_p_l = K * (density * velocity**2) / 2 # Local Pressure Drop Pa
delta_p = delta_p_f + delta_p_g + delta_p_l # Total Pressure Drop Pa

# Reynolds Number
Re = density * velocity * D_hydraulic / viscosity # Reynolds Number

# Mass Flow Rate
m_flow = density * velocity * A_flow # Mass Flow Rate kg/s

# Hydraulic Diameter
D_hydraulic = 4 * A_flow / P_wetted # Hydraulic Diameter [m]

# Prandtl Number c_p * mu / k
Pr = c_p * viscosity / th_conductivity # Prandtl Number

# Gas Conductivity
A_He = 15.8 # For Helium
A_Ar = 1.97 # For Argon
A_Kr = 1.15 # For Krypton
A_Xe = 0.72 # For Xenon
K_gas = A * 1e-4 * (T + 273.15) ** 0.79 # Gas Conductivity [W/mK]

# HTC in the gap
htc_conduction = K_gas / effective_thickness # Heat Transfer Coefficient in the gap [W/m2K]
htc_radiative = (4 * boltzmann * (T + 273.15) ** 3) / (1 / emissivity_1 + 1 / emissivity_2 - 1) # Radiative Heat Transfer Coefficient [W/m2K]
htc_contact = C * (conductivity_1 * conductivity_2) / (conductivity_1 + conductivity_2) * contact_pressure / (Hardness_Softer_material * sqrt(effective_thickness)) # Contact Heat Transfer Coefficient [W/m2K]
htc_gap = htc_conduction + htc_radiative + htc_contact # Total Heat Transfer Coefficient in the gap [W/m2K]


# Temp profile in the cladding RADIAL
T_cladding(r) = T_cladding_out - q_i / (2 * math.pi * k) * math.log(r / R_cladding_outer) # Temperature Profile in the cladding [K]

# Delta T inlet vs outlet coolant
T(z) = T_inlet + q_i_max / (m_flow * c_p) * sin(pi * z / extrapolated_length) * number_fuel_rods_per_channel # Delta T inlet vs outlet coolant [K]

# Heat generation in the fuel AXIAL
q_i = q_i_max * cos(pi * z / extrapolated_length) # Heat Generation in the fuel [W/m3]

extrapolated_length = length + 0.714 / Sigma_Transport # Extrapolated Length [m]

# Temp profile in the coolant RADIAL
T_coolant(r) = T_coolant_infinity + q_i / (2 * math.pi * k * r) # Temperature Profile in the coolant [K]

# Temp profile in the GAP RADIAL
T_gap(r) = T_gap_outer + q_i / (2 * math.pi * r * htc_gap) # Temperature Profile in the gap [K]

# Temp progile in the fuel RADIAL
void_factor = 1 - (ln((R_fuel_outer/r)**2))/((R_fuel_outer/r)**2 - 1)
K_fuel = correlation depending on T !
T_fuel(r) = T_fuel_outer + q_i / (4 * math.pi * k_fuel) * void_factor # Temperature Profile in the fuel [K]