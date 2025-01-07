import numpy as np
import nuclei_func as nu
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# properties correlations
def cp_cool(T):
    # coolant isobaric specific heat [J kg^-1 K^-1]
    T = T + 273.15  # K
    return 1608 - 0.7481 * T + 3.929e-4 * T**2

def mu_cool(T):
    # coolant dynamic viscosity [Pa s]
    T = T + 273.15  # K
    return np.exp(813.9/T - 2.530)*1e-3

def tc_cool(T):
    # coolant thermal conductivity [W m^-1 K^-1]
    T = T + 273.15  # K
    return 110 - 0.0648 * T + 1.16e-5 * T**2

def rho_cool(T):
    # coolant density [kg m^-3]
    T = (T * (9/5)) + 32  # °F
    return 954.1579 + T * (T * (T * 0.9667 * 1e-9 - 0.46 * 1e-5) - 0.1273534)


def tc_clad(T):
    # cladding thermal conductivity [W m^-1 K^-1]
    return 13.95 + 0.01163 * T


def tc_gap(T, x_He = 1):
    # gap thermal conductivity [W m^-1 K^-1]
    # x_He: He molar fraction [-], rest is considered Xe
    T = T + 273.15  # K
    return ( 15.8e-4 * T**0.79 )**x_He * ( 0.72e-4 * T**0.79 )**(1 - x_He)
    
def tc_fuel(T, x, y_Pu, p, Bu):
    # fuel thermal conductivity [W m^-1 K^-1]
    # T: temperature [°C]
    # x: O/M ratio [-]
    # y_Pu: Pu weight percent [-]
    # p: porosity [-]
    # Bu: burnup [GWd/tHM]

    T = T + 273.15  # K
    A = 0.01926 + 1.06e-6 * x + 2.63e-8 * y_Pu
    B = 2.39e-4 + 1.37e-13 * y_Pu
    D = 5.27e9
    E = 17109.5

    k0 = (1/(A + B*T) + D * np.exp(-E/T) / T**2)*(1 - p)**2.5

    k = 1.755 + (k0 - 1.755) * np.exp(-Bu/128.75)
    return k


def tc_fuel_integral(Ti, To, x, y_Pu, p, Bu):
    # fuel thermal conductivity integral [W m^-1]
    # Ti, To: boundary temperatures [°C]
    # x: O/M ratio [-]
    # y_Pu: Pu weight percent [-]
    # p: porosity [-]
    # Bu: burnup [GWd/tHM]
    # Ti, To = max(Ti, To), min(Ti, To)
    Ti = Ti + 273.15  #  K
    To = To + 273.15  # K
    A = 0.01926 + 1.06e-6 * x + 2.63e-8 * y_Pu
    B = 2.39e-4 + 1.37e-13 * y_Pu
    D = 5.27e9
    E = 17109.5

    ko_int = (np.log((A + B*Ti)/(A + B*To))/B + D *
                (np.exp(-E/Ti) - np.exp(-E/To))/E)*(1-p)**2.5

    k_int = 1.755*(Ti-To) + (ko_int - 1.755*(Ti-To)) * np.exp(-Bu/128.75)
    return k_int

def S_yield(T):
    """
    Calcola il valore di S_y,0.2% in base alla temperatura T (in °C).

    Args:
        T (float): Temperatura in gradi Celsius.

    Returns:
        float: Valore di S_y,0.2% (in MPa).
    """
    if T < 600:
        return 555.5 - 0.25 * T
    elif 600 <= T < 1000:
        return 405.5 - 0.775 * (T - 600)
    else:  # T >= 1000
        return 345.5 - 0.25 * T


def fission_rate(flux, rho_f):
    # computes the fission rate
    # flux: neutron flux [s^-1]
    # rho_f: fuel density [kg m^-3]
    # returns:
    # fission rate [s^-1]

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
    return Sigma_fuel * flux

def radial_redstr(R, r0):
    # computes the radial redistribution of Pu and power density
    # R: fuel outer radius [m]
    # r0: inner void radius [m]
    # returns:
    # chi: radial redistribution function
    # I: first radial moment integral function of chi
    
    D = 0.01
    alpha = 10
    J1 = R**2/(alpha) * ( (r0/R)*np.exp(-alpha*r0/R) - np.exp(-alpha) ) + R**2/(alpha**2) * ( np.exp(-alpha*r0/R) - np.exp(-alpha) )
    J2 = R**2/(2*alpha) * ( (r0/R)*np.exp(-2*alpha*r0/R) - np.exp(-2*alpha) ) + R**2/(4*alpha**2) * ( np.exp(-2*alpha*r0/R) - np.exp(-2*alpha) )

    r_star = R/alpha * np.log(2*J1/J2)
    chi = lambda r: 1 + D*( np.exp(-2*alpha*(r-r_star)/R) - 2*np.exp(-alpha*(r-r_star)/R) )
    I = lambda r: r**2/2 - D*R/(2*alpha) * r * np.exp(-2*alpha*(r - r_star)/R) - D*R**2/(4*alpha**2) * np.exp(-2*alpha*(r - r_star)/R) + 2*D*R/alpha * r * np.exp(-alpha*(r - r_star)/R) + 2*D*R**2/(alpha**2) * np.exp(-alpha*(r - r_star)/R)

    return chi, I


def thermal_up_to_fuel(q_lin, Rco, Rci, Rf, x_He = 1):
    # thermal analysis up to fuel outer radius
    # q_lin(z): linear power density [W m^-1]
    # Rco: cladding outer radius [m]
    # Rci: cladding inner radius [m]
    # Rf: fuel outer radius [m]
    # x_He: He molar fraction in the gap [-]
    # returns:
    # T_cool: coolant temperature [°C]
    # T_co: cladding outer temperature [°C]
    # T_ci: cladding inner temperature [°C]
    # T_gap: gap temperature [°C]
    # T_fo: fuel outer temperature [°C]

    # channel geometry
    f = 0.5
    pitch = 8.275e-3  # m
    Ap = 0.25*np.sqrt(3)*pitch**2 - 0.5*np.pi*Rco**2
    Dh = 4*Ap / Rco
    Ha = 0.85
    # number of axial points
    nz = len(Rco)
    z = np.linspace(-Ha/2, Ha/2, nz)

    # coolant
    G = 0.049 # kg s^-1
    def dTcool_dz(z, T):
        return f * q_lin(z) / (G * cp_cool(T))

    solution = solve_ivp(dTcool_dz, (-Ha/2, Ha/2), y0=395*np.ones(nz), t_eval=z)
    T_cool = solution.y[0]

    # convective HTC
    u = G / (Ap * rho_cool(T_cool))
    # u_max = np.max(u)
    Re = (rho_cool(T_cool) * u * Dh) / mu_cool(T_cool)
    Pr = (cp_cool(T_cool) * mu_cool(T_cool)) / tc_cool(T_cool)
    Pe = Re * Pr  # Pèclet number
    Nu = 7 + 0.025 * (Pe**(0.8))
    h = Nu * tc_cool(T_cool) / Dh

    # cladding outer
    T_co =  T_cool + f * q_lin(z) / (2 * np.pi * Rco * h)

    # cladding inner
    def fun_clad(T_ci):
        return T_ci - T_co - q_lin(z) * np.log(Rco/Rci) / (2 * np.pi * tc_clad(0.5*(T_ci + T_co)))

    T_ci = fsolve(fun_clad, T_co)
    # T_clad_fun = lambda r: (T_co - T_ci) * (r - Rci)/(Rco - Rci) + T_ci
    # r_cladding = np.linspace(Rci, Rco, 100)
    # T_clad = T_clad_fun(r_cladding)
    # Generare un array 2D per r_cladding
    r_cladding = np.array([np.linspace(Rci[i], Rco[i], 100) for i in range(len(T_ci))])
    T_clad = (T_co[:, None] - T_ci[:, None]) * (r_cladding - Rci[:, None]) / (Rco[:, None] - Rci[:, None]) + T_ci[:, None]

    # gap
    t_gap = Rci - Rf
    def fun_gap(T_fo):
        return T_fo - T_ci - q_lin(z) * t_gap / (2 * np.pi *Rf * tc_gap(0.5*(T_fo + T_ci), x_He))

    T_fo = fsolve(fun_gap, T_ci)
    T_gap = 0.5*(T_fo + T_ci)

    return T_cool, T_clad, T_co, T_ci, T_gap, T_fo

def thermal_in_fuel(q_lin, BU, Rf, T_fo, nr, restructured=False):
    # thermal analysis in the fuel
    # q_lin(z): linear power density [W m^-1]
    # BU: burnup [GWd tHM^-1]
    # Rf: fuel outer radius [m]
    # T_fo: fuel outer temperature [°C]
    # nr: number of radial nodes
    # restructured: restructured fuel flag
    # returns:
    # Tfuel: fuel temperature matrix
    # Rfuel: fuel radius matrix

    # active fuel height [m]
    Ha = 0.85
    # fuel theoretical density [kg m^-3]
    TD = 11310
    # fuel density [kg m^-3]
    rho_f = TD * 0.945 
    # initial Pu mass fraction
    y0 = 0.29
    # deviation from stoichiometry
    x = 2 - 1.945

    Rf_avg = np.mean(Rf)
    nz = len(Rf)
    z = np.linspace(-Ha/2, Ha/2, nz)

    if not restructured:
        # fuel cross section [m^2]
        Omega = np.pi*Rf_avg**2
        # fuel temperature matrix
        Tfuel = np.zeros((nz, nr))
        Rfuel = np.zeros((nz, nr))
        # thermal conductivity [W m^-1K^-1]
        kf = lambda T: tc_fuel(T, x, y0, 1 - 0.945, BU)
        for i in range(nz):
            r_fuel = np.linspace(Rf[i], 0, nr)

            def dTf_dr(r, T):
                return - q_lin(z[i]) * r /(2 * Omega* kf(T))

            solution = solve_ivp(
                dTf_dr, 
                (Rf[i], 0), 
                y0=[T_fo[i]], 
                t_eval=r_fuel
            )
            Tfuel[i, :] = np.flip(solution.y[0])
            Rfuel[i, :] = np.flip(solution.t)
    else:
        T1 = 1800
        T2 = 1600
        rho1 = 0.98 * TD
        rho2 = 0.96 * TD
        rho3 = rho_f
        p1 = 1 - rho1 / TD
        p2 = 1 - rho2 / TD
        p3 = 1 - rho3 / TD

        # piecewise thermal conductivity
        def kf(T):
            if T > T1:
                return tc_fuel(T, x, y0, p1, BU)
            elif T > T2:
                return tc_fuel(T, x, y0, p2, BU)
            else:
                return tc_fuel(T, x, y0, p3, BU)

        Tfuel = np.zeros((nz, nr))
        Rfuel = np.zeros((nz, nr))
        Tf_r = np.zeros(nr)
        r_f = np.zeros(nr)
        solution = []
        for i in range(nz):   
            # guess internal radius [m]
            r0 = 0.1 * Rf[i]

            toll = 1e-9
            max_iter = 1000
            err = toll + 1
            n = 0
            while err > toll and n <= max_iter:
                # fuel cross section [m^2]
                Omega = np.pi * (Rf[i]**2 - r0**2)
                # radial redistribution 
                chi, I = radial_redstr(Rf[i], r0)
                # piecewise power density along z [W m^-3]
                def q0(z, T):
                    if T > T1:
                        return (rho1/rho3) * q_lin(z) / Omega
                    elif T > T2:
                        return (rho2/rho3) * q_lin(z) / Omega
                    else:
                        return q_lin(z) / Omega
                
                def dTf_dr(r, T):
                    return - q0(z[i], T) / ( r * kf(T) ) * ( I(r) - I(r0) )

                sol = solve_ivp(dTf_dr, (Rf[i], r0), y0=[T_fo[i]], t_eval=np.linspace(Rf[i], r0, nr))
                Tf_r = np.flip(sol.y[0])
                r_f = np.flip(sol.t)

                # find r1, r2 where T = T1, T2
                idx = np.argmin(abs(Tf_r - T1))
                r1 = r_f[idx]
                idx = np.argmin(abs(Tf_r - T2))
                r2 = r_f[idx]

                # new void radius
                r0_new = np.sqrt( (rho1-rho2)*r1**2 / rho1 + (rho2-rho3)*r2**2 / rho1 )

                err = np.abs(r0_new - r0)
                r0 = r0_new
                n += 1

            # save values of r0, r1, r2, T0
            solution.append([r0, r1, r2])
            # update fuel temp matrix 
            Tfuel[i, :] = Tf_r
            Rfuel[i, :] = r_f

        r0 = np.array([sol[0] for sol in solution])
        r1 = np.array([sol[1] for sol in solution])
        r2 = np.array([sol[2] for sol in solution])
        print("Max radii and internal temperature")
        print(f"r0: {np.max(r0)*1000:.2f} mm")
        print(f"r1: {np.max(r1)*1000:.2f} mm")
        print(f"r2: {np.max(r2)*1000:.2f} mm")
        print(f"T0: {np.max(Tfuel[:,0]):.2f} °C")

        plt.figure()
        plt.title("Fuel restructuring radii")
        plt.plot(r0*1000, z, label='r0')
        plt.plot(r1*1000, z, label='r1')
        plt.plot(r2*1000, z, label='r2')
        plt.plot(Rf*1000, z, label='Rf')
        plt.xlabel('r [mm]')
        plt.ylabel('z [mm]')
        plt.legend()
        plt.grid()

    return Tfuel, Rfuel







# plotting functions

def plot_axial_T(T, r, r0, R, n_r, Ha):
    if r > R or r < r0:
        print("r not in range in plot_axial_T")
        return

    # Find the radial index corresponding to the given r
    radial_positions = np.linspace(r0, R, n_r)
    radial_index = np.argmin(abs(r - radial_positions))

    # Define the axial positions
    z = np.linspace(-Ha / 2, Ha / 2, T.shape[0])
    
    # Plot the temperature profile along the axial direction
    plt.figure()
    plt.plot(T[:, radial_index], z * 1000)
    plt.title(f"Fuel axial T [°C] @ r = {r * 100} mm")
    plt.xlabel("T [°C]")
    plt.ylabel("z [mm]")
    plt.grid()

def plot_radial_T(T, z, r0, R, n_r, Ha):
    # Find the axial index corresponding to the given z
    axial_positions = np.linspace(-Ha / 2, Ha / 2, T.shape[0])
    axial_index = np.argmin(abs(z - axial_positions))

    # Define the radial positions
    r = np.linspace(r0, R, n_r)

    # Plot the temperature profile along the radial direction
    plt.figure()
    plt.plot(r * 1000, T[axial_index, :])
    plt.title(f"Fuel radial T [°C] @ z = {z * 1000} mm")
    plt.xlabel("r [mm]")
    plt.ylabel("T [°C]")
    plt.grid()

def plot_Tfuel(Tfuel, Rfuel, Ha, restructured=False):
    
    
    # Create meshgrid for plotting
    z = np.linspace(-Ha/2, Ha/2, Tfuel.shape[0])
    r0_avg = np.mean(Rfuel[:,0])
    Rf_avg = np.mean(Rfuel[:,-1])
    r = np.linspace(r0_avg, Rf_avg, Tfuel.shape[1])
    R, Z = np.meshgrid(r, z, indexing='xy')
    
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Plot temperature colormap
    plt.pcolormesh(R*1000, Z*1000, Tfuel, shading='auto', cmap='hot')
    plt.colorbar(label='Temperature [°C]')
    
    if restructured:
        # Add contours for T=1800K and T=1600K
        CS1 = plt.contour(R*1000, Z*1000, Tfuel, levels=[1800], colors='white', linestyles='solid', linewidths=2)
        CS2 = plt.contour(R*1000, Z*1000, Tfuel, levels=[1600], colors='white', linestyles='solid', linewidths=2)
    
        # Add contour labels
        plt.clabel(CS1, inline=True, fmt='%1.0f°C')
        plt.clabel(CS2, inline=True, fmt='%1.0f°C')
    
    # Set labels and title
    plt.title('Fuel Temperature Distribution')
    plt.xlabel('r [mm]')
    plt.ylabel('z [mm]')
    plt.grid(True)
    plt.show()

def plot_T_colormap(T, r0, R, n_r, Ha):
    # Define the axial and radial positions
    z = np.linspace(-Ha / 2, Ha / 2, T.shape[0])
    r = np.linspace(r0, R, n_r)
    r, z = np.meshgrid(r, z, indexing="xy")

    # Plot the temperature colormap
    plt.figure()
    plt.title("Fuel temperature")
    plt.pcolormesh(r * 1000, z * 1000, T, shading='auto', cmap='hot')
    plt.colorbar(label='T [°C]')
    plt.xlabel('r [mm]')
    plt.ylabel('z [mm]')
    plt.grid()
    plt.show()




def plot_expanded_radii(z, Rf_hot, Rci_hot, Rco_hot, Rf_cold, Rci_cold, Rco_cold):
    plt.figure(figsize=(10, 6))
    plt.plot(Rf_hot * 1000, z, label='Rf_hot (mm)')
    plt.plot(Rci_hot * 1000, z, label='Rci_hot (mm)')
    plt.plot(Rco_hot * 1000, z, label='Rco_hot (mm)')
    plt.ylabel('Axial Position z (m)')
    plt.xlabel('Radii (mm)')
    plt.title('Expanded Radii Along Axial Position')
    plt.legend()
    plt.grid(True)
    plt.axvline(x=Rf_cold[0] * 1000, color='blue', linestyle='--', label='Rf_cold (mm)')
    plt.axvline(x=Rci_cold[0] * 1000, color='green', linestyle='--', label='Rci_cold (mm)')
    plt.axvline(x=Rco_cold[0] * 1000, color='red', linestyle='--', label='Rco_cold (mm)')
    plt.legend()
    plt.show()

def plot_axial_temperatures(z, Tcool, Tco, Tci, Tgap, Tf):
    plt.figure(figsize=(10, 6))
    plt.plot(Tcool, z, label='Coolant')
    plt.plot(Tco, z, label='Cladding Outer')
    plt.plot(Tci, z, label='Cladding Inner')
    plt.plot(Tgap, z, label='Gap')
    plt.plot(Tf[:, -1], z, label='fuel outer')
    plt.plot(Tf[:, 0], z, label='fuel inner')
    plt.xlabel('Axial Position z (m)')
    plt.ylabel('Temperature (°C)')
    plt.title('Axial temperature profiles')
    plt.legend()
    plt.grid(True)
    plt.show()
