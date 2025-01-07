# Fission Gas Release function to find the hignt of the planum:
import numpy as np
import nuclei_func as nu

def FGR(H_p, rho_f, Flux, Tf_max, Tg_max, t_irr, R_fi, R_fo, t_gap, H_a, p_lim):

    # We report the usefull data to find the concentrations, 
    R_grain = 1e-5 #[m]
    y_Fg = 0.30
    sigmaU235_f = 1.047756375 #[b]
    sigmaU238_f = 0.55801001 #[b]
    sigmPu_f = 1.689844625 #[b]
    S_U235 = nu.macro(sigmaU235_f, rho_f, 235)
    S_U238 = nu.macro(sigmaU238_f, rho_f, 238)
    S_Pu = nu.macro(sigmPu_f, rho_f, 239)
    quantity = [S_U235, S_U238, 0, S_Pu] 
    mol_quality = nu.w2mol([0.711, 0.29], [235 + 32, 239 + 32])
    quality = [mol_quality[0]*(1 - mol_quality[1]), (1 - mol_quality[0])*(1 - mol_quality[1]), 2, mol_quality[1]]
    Sigma_fuel = nu.mixture(quantity, quality, 1)

    # And we usef the found cross sections to find the fission rate (F)
    F = Sigma_fuel * Flux

    # We define the effective Diffusion coefficient (D):
    Tf_max_K = Tf_max + 273.15 #[K]
    D = 5 * 1e-8 * np.exp(- 40262 / Tf_max_K)
    DDD = 5 * 1e-8 * np.exp(- 40262 /1673)
    print(f"Diffusion coefficient: {DDD}")

    # We have solved the time dependent diffusion equation for the fission gasses to obtain:
    V_grain = (4 * np.pi * R_grain**3) / 3
    p = y_Fg * F * t_irr * V_grain
    g = 8 * y_Fg * F * R_grain**5/(np.pi**3 * D) * (1 - np.exp(- (D * t_irr) * (np.pi / R_grain)**2))
    print(f"Release fraction: {1 - g/p:.9f}")
    print(f"g: {g:.9f}, p: {p:.9f}")
    R_Fg = p - g
    # R_Fg = y_Fg * F * (t_irr*V_grain - ((8 * y_Fg * (R_grain**5))/((np.pi**3) * D)) * (1 - np.exp(- (D * t_irr) * (np.pi / R_grain)**2)))
    V_f = np.pi * (R_fo**2 - R_fi**2) * H_a
    R_tot = R_Fg * (V_f / V_grain) # total nuber of Fg atoms released in the plenum

    # Now we caluclate the plenum size to accomodate all Fg's:
    p_He = 0.1 * 1e6 #[Pa]
    T_He = 20 + 273.15 #[K]
    V_p = np.pi * (R_fo + t_gap)**2 * H_p #[m^3]
    n_He = p_He * V_p / (8.31446 * T_He)

    # We add the fission gasses to the plenum and find the final pressure
    n_Fg = R_tot / (6.022 * 1e23)
    n_tot = n_He + n_Fg
    Tg_max_K = Tg_max + 273.15 #[K]
    p_tot = (n_tot * 8.31446 * Tg_max_K / V_p) * 1e-6 #[MPa]

    # Finally we compare the obtained pressure to the limit
    if p_tot < p_lim:
        print('Hell yeah that plenum is on point!!!')
    else:
        print('What the actual sigma... ew')

    # We also find the molar fraction for the k_gap evluation 
    y_Fg = n_Fg / n_tot
    A_Fg = np.mean((1.15, 0.72))

    return y_Fg, A_Fg, p_tot, n_Fg, n_He
