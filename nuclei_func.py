############################################################################################################
# NUCLEI FUNCTIONS
#
# This file contains a collection of functions related to nuclear physics and reactions.
# Each function serves a specific purpose and can be used independently or in combination with others.
# The functions in this file are designed to perform calculations and provide useful information
# for studying and analyzing nuclear phenomena.
#
# Author: [Pagliuca Simone]
# Date: [18/04/2024]
#
############################################################################################################
import inspect
import numpy as np


############################################################################################################
# Function to calculate the macroscopic cross section from the microscopic cross section
############################################################################################################
def macro(micro_sigma, density, molar_mass):
    N_a = 6.022e23 # Avogadro's number
    barn = 1e-24 # b-->cm2
    macro_sigma = micro_sigma * N_a * density / molar_mass * barn
    return macro_sigma

############################################################################################################
# Function to calculate the cross section of a mixture
############################################################################################################
def mixture(quantity, quality, normalization_cond='no_normalize'):

    #check if the quantity and quality lists have the same length
    if len(quantity) != len(quality): 
        print("The quantity and quality lists must have the same length.")
        return None

    #check if the sum of the quality values is equal to 1, if not, ask the user if they want to normalize them
    if sum(quality) != 1:
        if normalization_cond == 'normalize':
            quality = [q/sum(quality) for q in quality]
        elif normalization_cond == 'no_normalize':
            print(f"\033[93mThe sum of the qualites values is not 1, at line {inspect.currentframe().f_back.f_lineno}, default behaviour is ignoring.\033[0m")  
        elif normalization_cond == 'silent':
            pass
    #calculate the macroscopic cross section of the mixture
    temp = [quantity[i]*quality[i] for i in range(len(quantity))]
    mixture_val = sum(temp)
    return mixture_val

############################################################################################################
# From weight percentage to molar fraction
############################################################################################################
def w2mol(weight_fraction, molar_mass):

    # check if the weight percentage and atomic weight lists have the same length
    if len(weight_fraction) != len(molar_mass):
        print("The weight percentage and atomic weight lists must have the same length.")
        return None
    
    ans = [] #initialize the molar fraction list
    denominator = sum([weight_fraction[jj]/molar_mass[jj] for jj in range(len(weight_fraction)) ]) #compute the denominator outside the loop
    #loop to calculate the molar fractions
    for ii in range(len(molar_mass)):
        ans.append((weight_fraction[ii]/molar_mass[ii])/denominator)

    return ans

############################################################################################################
# From molar fraction to weight percentage
############################################################################################################
def mol2w(mol_fraction, molar_mass):

    # check if the weight percentage and atomic weight lists have the same length
    if len(mol_fraction) != len(molar_mass):
        print("The weight percentage and atomic weight lists must have the same length.")
        return None
    
    ans = [] #initialize the weight fraction list
    denominator = sum([mol_fraction[jj]*molar_mass[jj] for jj in range(len(mol_fraction)) ]) #compute the denominator outside the loop
    #loop to calculate the molar fractions
    for ii in range(len(molar_mass)):
        ans.append((mol_fraction[ii]*molar_mass[ii])/denominator)

    return ans

############################################################################################################
# From volumetric percentage to weight fraction
############################################################################################################
def vol2w(volume_fraction, density):

    # check if the weight percentage and atomic weight lists have the same length
    if len(volume_fraction) != len(density):
        print("The weight percentage and atomic weight lists must have the same length.")
        return None
    
    ans = [] #initialize the molar fraction list
    denominator = sum([volume_fraction[jj]*density[jj] for jj in range(len(volume_fraction))]) #compute the denominator outside the loop
    #loop to calculate the molar fractions
    for ii in range(len(density)):
        ans.append((volume_fraction[ii]*density[ii])/denominator)

    return ans

############################################################################################################
# This is used for exercises 7 and 8
############################################################################################################
def compute_k(materials):
    densities = [ii['density'] for ii in materials]
    # print("densities",densities)
    volume_fractions = [ii['vol_fraction']*(0.3*ii['fuel'] + 1*(not ii['fuel'])) for ii in materials]
    # print("volume fract", volume_fractions)
    weight_fractions = vol2w(volume_fractions, densities)
    # print("weight fractions", weight_fractions)
    weight_fractions_fuel = [weight_fractions[ii]*(materials[ii]['fuel']) for ii in np.arange(len(materials))]
    # print("weight fraction fuel", weight_fractions_fuel)
    Macro_abs = [macro(ii['sigma_a'], ii['density'], ii['molar_mass']) for ii in materials]
    print("Macro_abs", Macro_abs)
    Macro_Fiss = [macro(ii['sigma_f'], ii['density'], ii['molar_mass']) for ii in materials]
    print("Macro_Fiss", Macro_Fiss)
    Nu_Sigma_fiss = [materials[ii]['nu'] * Macro_Fiss[ii] for ii in np.arange(len(materials))]
    # print("Nu_Sigma_fiss", Nu_Sigma_fiss)

    ############################################################################################################
    # Multiplication factor K_inf
    ############################################################################################################
    eta = mixture(Nu_Sigma_fiss, weight_fractions, 'silent') / mixture(Macro_abs, weight_fractions_fuel, 'silent')
    f = mixture(Macro_abs, weight_fractions_fuel, 'silent') / mixture(Macro_abs, weight_fractions, 'silent')
    p = 1 #non-leakage probability
    epsilon = 1

    K_inf = eta * epsilon * p * f
    A_mass_percent = round(weight_fractions_fuel[0] * 100, 2)
    B_mass_percent = round(weight_fractions_fuel[1] * 100, 2)
    return K_inf,A_mass_percent, B_mass_percent

############################################################################################################
# This is used for exercise 6
############################################################################################################
def six_factors(compounds, densities, materials, volumes):
    # Initialize lists
    Macro_Fission = []
    Macro_Absorption = []
    Macro_Fission_Nu = []

    # Work on each compound, one by one
    for ii, compound in enumerate(compounds):
        # Initialize lists
        sigma_f_list = []
        sigma_a_list = []
        atom_fractions = []
        molar_masses = []
        nu_list = []

        # Get data for each component in the compound from the materials list
        for component in compound:
            for material in materials:
                if material['name'] == component['name']:
                    sigma_f_list.append(material['sigma_f'])
                    sigma_a_list.append(material['sigma_c']+material['sigma_f']) # Absorption = Capture + Fission
                    nu_list.append(material['nu'])
                    atom_fractions.append(component['atom_fraction'])
                    molar_masses.append(material['molar_mass'])
                    break
        # Macroscopic cross sections for each element in the compound, using compound density
        macroscopic_sigma_f = [macro(sigma_f_list[jj], densities[ii], molar_masses[jj]) for jj in np.arange(len(sigma_f_list))]
        macroscopic_sigma_a = [macro(sigma_a_list[jj], densities[ii], molar_masses[jj]) for jj in np.arange(len(sigma_a_list))]
        nu_macro_f = [nu_list[jj] * macroscopic_sigma_f[jj] for jj in np.arange(len(nu_list))]
        
        # Overall Macroscopic cross sections for the compound
        macroscopic_sigma_f = mixture(macroscopic_sigma_f, atom_fractions, 'silent')
        macroscopic_sigma_a = mixture(macroscopic_sigma_a, atom_fractions, 'silent')
        nu_f = mixture(nu_macro_f, atom_fractions, 'silent')

        # Append the values to the lists
        Macro_Fission.append(macroscopic_sigma_f)
        Macro_Absorption.append(macroscopic_sigma_a)
        Macro_Fission_Nu.append(nu_f)

        print("\n")

    ############################################################################################################
    # Compute the parameters for the 6 factor formula
    ############################################################################################################
    # When we have to use fractions for the fuel only we get rid of the last value
    # This could be done in a fancier way by adding a flag to each component 
    # but since it's only one we are going to use
    # as a convention that the last is the moderator or coolant

    # Compute fraction of each compound in the core
    tot_vol = sum(volumes)
    vol_fractions = [volume / tot_vol for volume in volumes]
    print("vol_fractions", vol_fractions)

    # Convert into weight fractions
    fractions = vol_fractions #vol2w(vol_fractions, densities)

    eta = sum([Macro_Fission_Nu[ii]/Macro_Absorption[ii] for ii in np.arange(len(Macro_Fission_Nu)-1)])
    f = mixture(Macro_Absorption[:-1], fractions[:-1], 'silent') / mixture(Macro_Absorption, fractions, 'silent')
    return eta, f