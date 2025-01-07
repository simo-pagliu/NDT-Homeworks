import numpy as np 
import functions as fn 

def main(t_gap, H_plenum): 
    # fast neutron flux ( E > 100 keV ) [n cm^-2 s^-1]
    neutron_flux = 6.1e15

    # number of axial nodes
    nz = 25
    # number of radial nodes in the fuel
    nr = 100

    # active height [m]
    Ha = 0.85

    # cold geometry
    Rf = 5.42e-3/2 * np.ones(nz)
    Rci = Rf + t_gap
    Rco = 6.55e-3/2 * np.ones(nz)

    V_fuel = np.pi*Ha*Rf[0]**2
    TD = 11310
    rho_f = 0.945 * TD

    # irradiation time [d]
    t_irr = 365

    # axial linear power evaluation
    q_0 = 38.7e3  # W m^-1
    c = np.array([0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.49])
    z_fit = np.array([42.5, 127.5, 212.5, 297.5, 382.5, 467.5,552.5, 637.5, 722.5, 807.5])*0.001 - Ha/2 # m
    q_fit = q_0 * c
    polyn = np.poly1d(np.polyfit(z_fit, q_fit, q_fit.size-1))
    z = np.linspace(-Ha/2, +Ha/2, nz)
    def q_lin(z): return polyn(z)  # q = q(z) [W m^-1]
    # total rod power [W]
    Q = np.trapz(q_lin(z), z)
    # burnup [GWd tHM^-1]
    BU = 1e-6 * t_irr * Q / (V_fuel * rho_f)

    print(f"Total rod power: {Q:.3f} W")
    print(f"Burnup at {t_irr} d: {BU:.3f} GWD/tHM")


            # 1. Thermal analysis at cold geometry

    Tcool, Tco, Tci, Tgap, Tfo = fn.thermal_up_to_fuel(q_lin, Rco, Rci, Rf)
    Tfuel, Rfuel = fn.thermal_in_fuel(q_lin, BU, Rf, Tfo, nr, restructured=False)
    # fn.plot_axial_temperatures(z, Tcool, Tco, Tci, Tgap, Tfuel)

    print("\nCold geometry")
    print(f"Rf: {Rf[0]*1000} mm, Rci: {Rci[0]*1000} mm, Rco: {Rco[0]} mm")
    print("Max temperatures [°C]")
    print(f"fuel: {np.max(Tfuel):.2f}, gap: {np.max(Tgap):.2f}, cladding: {np.max(Tci):.2f}, coolant: {np.max(Tcool):.2f}")
    Tf_max = np.max(Tfuel)
            
            # 2. FGR and plenum dimensioning

    R_grain = 1e-5 #[m]
    V_grain = 4/3 * np.pi * R_grain**3
    y_Fg = 0.30

    # fission rate [# s^-1]
    F_rate = fn.fission_rate(neutron_flux, rho_f)
    print(f"\nFission rate: {F_rate:.3e}")

    # fission gas production [atoms]
    p = lambda t: y_Fg * F_rate * t * V_grain

    # effective diffusion coefficient [m^2 s^-1]
    # conservative approach: max fuel temperature at cold geometry (hottest fuel)
    Tf_max_K = Tf_max + 273.15 #[K]
    D = 5 * 1e-8 * np.exp(- 40262 / Tf_max_K)
    # print(f"Deff: {D:.2e} m^2 s^-2")

    # fission gas retention in the grain [atoms]
    g = lambda t: y_Fg * F_rate * 8 * R_grain**5 / (np.pi**3 * D) * (1 - np.exp(-np.pi**2 * D * t / R_grain**2))

    # release, by all the grains
    FG_release = lambda t: ( p(t) - g(t) ) * V_fuel / V_grain

    # fractional release
    Frac_release = lambda t: ( p(t) - g(t) ) / p(t)
    print(f"Fractional release at {t_irr} d: {Frac_release(t_irr*24*3600):.10f}")

    # initial He moles
    p_He = 101325 #[Pa]
    T_He = 20 + 273.15 #[K]
    V_plenum = np.pi * Rci[0]**2 * H_plenum
    V_gas = V_plenum + np.pi * t_gap**2 * Ha
    n_He = p_He * V_gas / (8.31446 * T_He)
    # fission gas moles
    NA = 6.022e23
    print(f"FG release at {t_irr} d: {FG_release(t_irr):.3e}")
    n_fg = FG_release(t_irr) / NA
    x_He = n_He / (n_He + n_fg)
    print(f"He molar fraction: {x_He}")
    print(f"Thermal conductivity He only: {fn.tc_gap(np.max(Tgap))}")
    print(f"Thermal conductivity with FG: {fn.tc_gap(np.max(Tgap), x_He)}")

    # plenum pressure
    Tgap_max_K = np.max(Tgap) + 273.15
    p_tot = (n_He + n_fg) * 8.31446 * Tgap_max_K / V_plenum
    print(f"Plenum pressure: {p_tot/1e6:.3f} MPa")

if __name__ == '__main__':
    t_gap = 0.16e-3
    H_plenum = 0.3
    main(t_gap, H_plenum)
