# import of libraries

import numpy as np
import pandas as pd
from scipy import interpolate as intp
from scipy import integrate as intg
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

# parameters

global prm

prm = { # settings
        "series":True,                  # if multiple simulations are to be run in a row
        "graphs":True,                  # if graphs are to be plotted
        "long_graphs":False,             # if long term graph are to be plotted
        # simulation parameters
        "dt":0.0010,                    # [s] simulation time step (set to 0.0010 for standard precision, and to 0.0001 for increased precision)
        "t":2.0,                        # [s] simulation duration
        "f":2.0,                        # [Hz] liquid switch frequency
        # configuration parameters
        "Ngd":1,                        # [] number of gadoilinium bars in the generator, "1", "2" and "3" are allowed values, defaults to "1" if disallowed values are selected
        "h":0.01,                       # [m] distance between ferromagnetic plates                          
        "Aa":0.0016,                    # [m^2] surface area between plates exposed to air
        # thermodynamic parameters
        "Tc":290.0,                     # [K] cold liquid temperature
        "Th":300.0,                     # [K] hot liquid temperature
        "K":160.0,                      # [W/K] heat transfer parameter
        # gadolinium parameters
        "rho":7900.0,                   # [kg/m^3] gadolinium density
        "mu":0.15725,                   # [kg/mol] gadolinium molar mass
        "Ag":0.0004,                    # [m^2] surface area between plates exposed to one gadolinium bar
        # magnet parameters
        "Br":1.4,                       # [T] magnet remanence
        "Am":0.0004,                    # [m^2] surface area between plates exposed to magnet
        # coil parameters
        "Nc":10,                       # [] number of coil loops 
        "Rc":2000,                      # [Ohm] one coil resistance
        # load parameters
        "Rl":2000,                      # [Ohm] load resistance
        "Cl":1.0e+5,                    # [F] load capacitance
        # constants
        "mu0":1.257e-6                  # [N/A^2] magnetic constant
      }

if prm["Ngd"] not in [1, 2, 3]:
    prm["Ngd"] = 1

# series setup

# number of simulations:
n_series = 20
# changing parameters [name, start, end]:
series = [["Rl", 500, 4500]]
# x-axis label:
series_x = "Resistance ratio"
# x-axis range:
range_x = np.linspace(0.25, 2.25, n_series)

states_setup = []
linsp = []
for param in series:
    linsp.append(np.linspace(param[1], param[2], n_series))
for n in range(n_series):
    lst_par = [[], []]
    for m in range(len(series)):
        lst_par[0].append(series[m][0])
        lst_par[1].append(linsp[m][n])
    states_setup.append(lst_par)

# interpolation and testing of data

def read_data(path):
    df = pd.read_csv(path)
    prts = path.split("/")
    prts = prts[-1].split(",")
    name = prts[-1][:-4]
    try:
        val = float(prts[0])
    except:
        val = None
    df[name] = val
    return df

def interb(data):
    return intp.SmoothBivariateSpline(np.array(data)[:, 0], np.array(data)[:, 2], np.array(data)[:, 1])

def plot_rb(data, cc, start, stop, lstd, unitd):
    lstc = ['red', 'orange', 'yellow', 'lime', 'cyan', 'blue', 'purple']
    for i in range(len(lstd)):
        try:
            clr = lstc[i]
        except:
            clr = 'black'
        plt.plot(np.arange(start, stop, (stop-start)/300), 
                 cc*data(np.arange(start, stop, (stop-start)/300), lstd[i]), 
                 color=clr, label=f'${lstd[i]}{unitd}$')
    plt.grid()

def tb_test(rang):
    tb_calc = np.array([])
    f00 = lambda T : hc(T, 0.0)
    f75 = lambda T : hc(T, 2.0)
    F00 = lambda T : intg.quad(f00, 0, T)[0]
    F75 = lambda T : intg.quad(f75, 0, T)[0]
    for T00 in rang:
        x1 = 300.0
        n = 0
        eta = 1.0
        while eta >= 1.0e-6 and n < 1200:
            x0 = x1
            x1 = x0 - (F75(x0) - F00(T00))/(f75(x0))
            eta = abs(x1 - x0)
            n += 1
        tb_calc = np.append(tb_calc, (x1-T00))
    return tb_calc

def rainbow_test():
    plot_rb(hc, (prm["rho"]/prm["mu"]), 200, 350, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5], 'B')
    plt.xlabel('Temperature $[K]$')
    plt.ylabel('Heat capacity $\\big[\\frac{J}{m^{3} \\cdot K}\\big]$')
    plt.title('Gadolinium Heat Capacitance')
    plt.legend()
    plt.savefig("hc.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    plot_rb(m0, (prm["rho"]), 0, 5, [250, 260, 270, 280, 290, 300, 310], 'K')
    plt.xlabel('Magnetic induction $[T]$')
    plt.ylabel('Magnetization $\\big[\\frac{A}{m}\\big]$')
    plt.title('Gadolinium Magnetization Along $0001$ Axis')
    plt.legend()
    plt.savefig("m0.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    plot_rb(m1, (prm["rho"]), 0, 5, [250, 260, 270, 280, 290, 300, 310], 'K')
    plt.xlabel('Magnetic induction $[T]$')
    plt.ylabel('Magnetization $\\big[\\frac{A}{m}\\big]$')
    plt.title('Gadolinium Magnetization Along $1010$ Axis')
    plt.legend()
    plt.savefig("m1.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    plt.plot(np.arange(280, 320, 0.1), tb(np.arange(280, 320, 0.1)), color='#27a7d8', label='Original')
    plt.plot(np.arange(280, 320, 0.1), tb_test(np.arange(280, 320, 0.1)), color='#ff9a00', label='Recalculated')
    plt.grid()
    plt.xlabel('Initial temperature at $B=0.0T$ $[K]$')
    plt.ylabel('Temperature change when $B=2.0T$ $[K]$')
    plt.title('Gadolinium Magnetocaloric Effect')
    plt.legend()
    plt.savefig("tb.pdf", format="pdf", bbox_inches="tight")
    plt.show()

heat_capacity = pd.concat([read_data("Data/heat_capacity/0.0,T.csv"),
                           read_data("Data/heat_capacity/2.0,T.csv"),
                           read_data("Data/heat_capacity/5.0,T.csv"),
                           read_data("Data/heat_capacity/7.5,T.csv"),
                           read_data("Data/heat_capacity/10.0,T.csv")])

magnetization_0001 = pd.concat([read_data("Data/magnetization_0001/237.0,K.csv"),
                                read_data("Data/magnetization_0001/247.2,K.csv"),
                                read_data("Data/magnetization_0001/267.6,K.csv"),
                                read_data("Data/magnetization_0001/277.8,K.csv"),
                                read_data("Data/magnetization_0001/288.1,K.csv"),
                                read_data("Data/magnetization_0001/298.4,K.csv"),
                                read_data("Data/magnetization_0001/318.9,K.csv"),
                                read_data("Data/magnetization_0001/324.0,K.csv")])

magnetization_1010 = pd.concat([read_data("Data/magnetization_1010/236.9,K.csv"),
                                read_data("Data/magnetization_1010/247.1,K.csv"),
                                read_data("Data/magnetization_1010/267.6,K.csv"),
                                read_data("Data/magnetization_1010/277.4,K.csv"),
                                read_data("Data/magnetization_1010/288.0,K.csv"),
                                read_data("Data/magnetization_1010/298.4,K.csv"),
                                read_data("Data/magnetization_1010/318.8,K.csv"),
                                read_data("Data/magnetization_1010/324.0,K.csv")])

dt_at_db = read_data("Data/dt_at_db/2.0,T.csv")

global hc, m0, m1, tb

hc = interb(heat_capacity)              # x : temperature [K]; y : magnetic field [T]; z : heat capacity [J/mol K]
m0 = interb(magnetization_0001)         # x : magnetic field [T]; y : temperature [K]; z : magnetization [emu/g]
m1 = interb(magnetization_1010)         # x : magnetic field [T]; y : temperature [K]; z : magnetization [emu/g]

tb = intp.CubicSpline(np.array(dt_at_db)[:, 0], np.array(dt_at_db)[:, 1])            # x : temperature [K]; y : magnetocaloric effect at 7.5T pulsed magnetic field [K]


# functions


def gd_heat(Bg, Tg):
    '''
    Finds heat stored in gadolinium bar at specified magnetic induction and temperature.

    Parameters
    ----------
    Bg : float
        Magnetic induction in gadolinium bar, [T].
    Tg : float
        Temperature of gadolinium bar, [K].

    Returns
    -------
    Q : float
        Heat stored in gadolinium, [J].

    '''
    hcf = lambda T : prm["h"]*prm["Ag"]*(prm["rho"]/prm["mu"])*hc(T, Bg)[0][0]
    Q = intg.quad(hcf, 0, Tg)[0]
    return Q

def gd_temp(Bg, Q):
    '''
    Finds temperature of gadolinium bar at specified magnetic induction from stored heat.

    Parameters
    ----------
    Bg : float
        Magnetic induction in gadolinium bar, [T].
    Q : float
        Heat stored in gadolinium, [J].

    Returns
    -------
    Tg : float
        Temperature of gadolinium bar, [K].

    '''
    hcf = lambda T : prm["h"]*prm["Ag"]*(prm["rho"]/prm["mu"])*hc(T, Bg)[0][0]
    x1 = 300.0
    n = 0
    eta = 1.0
    while eta >= 1.0e-6 and n < 1200:
        x0 = x1
        x1 = x0 - (intg.quad(hcf, 0, x0)[0] - Q)/(prm["h"]*prm["Ag"]*(prm["rho"]/prm["mu"])*hc(x0, Bg)[0][0])
        eta = abs(x1 - x0)
        n += 1
    Tg = x1
    return Tg

def gd_temp_step(Bg0, Tg0, Bg1, Tl):
    '''
    Finds new temperature of a gadolinium bar.

    Parameters
    ----------
    Bg0 : float
        Initial magnetic induction in gadolinium bar, [T].
    Tg0 : float
        Initial temperature of gadolinium bar, [K].
    Bg1 : float
        Final magnetic induction in gadolinium bar, [T].
    Tl : float
        Temperature of liquid in heat exchanger, [K].

    Returns
    -------
    Tg1 : float
        Final temperature of gadolinium bar, [K].

    '''
    Q = gd_heat(Bg0, Tg0) + prm["dt"]*prm["K"]*(Tl - Tg0)
    Tg1 = gd_temp(Bg1, Q)
    return Tg1

def gd_magn(Bg, Tg):
    '''
    Finds averaged magnetization in gadolinium bar.

    Parameters
    ----------
    Bg : float
        Magnetic induction in gadolinium bar, [T].
    Tg : float
        Temperature of gadolinium bar, [K].

    Returns
    -------
    M : float
        Magnetization in gadolinium bar, [A/m].

    '''
    if Bg == 0:
        M = 0
    elif Bg > 0:
        M = ((prm["rho"])*m0(Bg, Tg)[0][0] + (prm["rho"])*m1(Bg, Tg)[0][0])/2
    elif Bg < 0:
        M = -((prm["rho"])*m0(-Bg, Tg)[0][0] + (prm["rho"])*m1(-Bg, Tg)[0][0])/2
    return M

def gd_b_field(Bg0, Tg, I):
    '''
    Finds total and electromotive magnetic fields in every gadolinium bar and corresponding coil.

    Parameters
    ----------
    Bg0 : list of floats
        List of initial total magnetic inductions in gadolinium bars, [T].
    Tg : list of floats
        List of gadolinium bar temperatures, [K].
    I : list of floats
        List of electric currents in every coil, [A].

    Returns
    -------
    Bg1 : numpy.array of floats
        List of final total magnetic inductions in gadolinium bars, [T].
    Bemf : numpy.array of floats
        List of final electromotive magnetic inductions in gadolinium bars, [T].

    '''
    Bg1 = np.zeros(prm["Ngd"])
    Bemf = np.zeros(prm["Ngd"])
    Bmag = np.zeros(prm["Ngd"])
    Bcur = np.zeros(prm["Ngd"])
    Bmut = np.zeros(prm["Ngd"])
    murmu0 = np.zeros(prm["Ngd"])
    for i in range(prm["Ngd"]):
        M = gd_magn(Bg0[i], Tg[i])
        if M == 0:
            murmu0[i] = prm["mu0"]
        else:
            murmu0[i] = abs(Bg0[i]/(Bg0[i]/prm["mu0"] - M))
            if murmu0[i] == 0.0 or murmu0[i] != murmu0[i]:
                murmu0[i] = prm["mu0"]
    Bmag[0] = prm["mu0"]*gd_magn(Bg0[0], Tg[0])*(prm["Aa"] + prm["Am"] + prm["Ag"]*(prm["Ngd"]-1)) + prm["Br"]*prm["Am"]
    for i in range(1, prm["Ngd"]):
        Bmag[0] -= prm["mu0"]*gd_magn(Bg0[i], Tg[i])*prm["Ag"]
    Bmag[0] /= (prm["Aa"] + prm["Am"] + prm["Ag"]*prm["Ngd"])
    for i in range(1, prm["Ngd"]):
        Bmag[i] = Bmag[0] - prm["mu0"]*gd_magn(Bg0[0], Tg[0]) + prm["mu0"]*gd_magn(Bg0[i], Tg[i])
    for i in range(prm["Ngd"]):
        Bcur[i] = (murmu0[i]*I[i]*prm["Nc"]**2)/(prm["h"])
        for j in range(prm["Ngd"]):
            if j == i:
                pass
            else:
                Bmut[i] += murmu0[i]*I[j]*prm["Nc"]/prm["h"] - murmu0[i]*Bcur[j]/murmu0[j]
        Bg1[i] = Bmag[i] + Bcur[i] + Bmut[i]
        Bemf[i] = Bmag[i] + Bmut[i]
    return Bg1, Bemf, Bmag

def coil_volt(Bg0, Bg1):
    '''
    Finds electromotive force in coil around gadolinium bar.

    Parameters
    ----------
    Bg0 : float
        Initial magnetic induction in gadolinium bar, [T].
    Bg1 : float
        Final magnetic induction in gadolinium bar, [T].

    Returns
    -------
    emf : float
        Electromotive force, [V].

    '''
    emf = - prm["Nc"]*(Bg1 - Bg0)*prm["Ag"]/(prm["dt"])
    return emf

def circuit_state(Tg, Bg, emf, I0, q0):
    '''
    Finds circuit state from gadolinium bar temperature, magnetic indiction in gadolinium bar, electromotive force and previous state.

    Parameters
    ----------
    Tg : float
        Temperature of gadolinium bar, [K].
    Bg : float
        Magnetic induction in gadolinium bar, [T].
    emf : float
        Electromotive force, [V].
    I0 : float
        Initial electric current, [A].
    q0 : float
        Initial capacitor charge, [F].

    Returns
    -------
    I1 : float
        Final electric current, [A].
    q1 : float
        Final capacitor charge, [F].

    '''
    M = gd_magn(Bg, Tg)
    if M == 0:
        murmu0 = prm["mu0"]
    else:
        murmu0 = Bg/(Bg/prm["mu0"] - M)
        if murmu0 == 0 or murmu0 != murmu0:
            murmu0 = prm["mu0"]
    L = abs((murmu0*prm["Ag"]*prm["Nc"]**2)/(prm["h"]))
    diff = lambda Y, t: [Y[1], (emf - (prm["Rc"] + prm["Rl"])*Y[1] - Y[0]/prm["Cl"])/L]
    q1s, I1s = intg.odeint(diff, [q0, I0], [0.0, prm["dt"]])[-1]
    return I1s, q1s

def gen_out(U, I):
    '''
    Finds output power from current and voltage.

    Parameters
    ----------
    U : float
        Total voltage from coil, [V].
    I : float
        Electric current, [A].

    Returns
    -------
    P : float
        Output power, [W].

    '''
    P = U*I*prm["Rl"]/(prm["Rl"]+prm["Rc"])
    return P

# series

if prm["series"]:
    states = states_setup
else:
    states = [[["series"], [False]]]

# main values

def main_values(t_data, T_data, U_data, I_data, q_data, P_data, Bt_data):
    indx = int(1/(prm["f"]*prm["dt"]))
    T_max = max(T_data[-indx:-1])
    T_min = min(T_data[-indx:-1])
    U_eff = 0
    for u in U_data[-indx:-1]:
        U_eff += u**2
    U_eff = (prm["Rl"]/(prm["Rl"]+prm["Rc"]))*np.sqrt(U_eff/indx)
    P_eff = sum(P_data[-indx:-1])/indx
    B_avg = np.mean(Bt_data)
    return [T_max, T_min, U_eff, P_eff, B_avg]

# simulation

n_sim = 1
results = []
for state in states:
    print(f'\n\n\nSimulation {n_sim} of {len(states)}')
    n_sim += 1
    for s in range(len(state[0])):
        prm[state[0][s]] = state[1][s]
    prm["Ngd"] = int(prm["Ngd"])
    t_data = np.arange(0, prm["t"], prm["dt"])              # time data
    T_data = np.zeros([len(t_data), prm["Ngd"]])            # temperature data
    U_data = np.zeros([len(t_data), prm["Ngd"]])            # coil electromotive force data
    I_data = np.zeros([len(t_data), prm["Ngd"]])            # output current data
    q_data = np.zeros([len(t_data), prm["Ngd"]])            # capacitor charge data
    P_data = np.zeros([len(t_data), prm["Ngd"]])            # output power data
    Bt_data = np.zeros([len(t_data), prm["Ngd"]])           # totalt magnetic induction data
    BE_data = np.zeros([len(t_data), prm["Ngd"]])           # electromotive magnetic induction data
    Bm_data = np.zeros([len(t_data), prm["Ngd"]])           # magnet magnetic induction data
    Tl = np.zeros(prm["Ngd"])
    t_switch = 1/prm["f"]/prm["Ngd"]/2
    n_switch = 0
    l_state = [False, True, False]
    for g in range(prm["Ngd"]):
        T_data[0, g] = (prm["Th"]+prm["Tc"])/2
        if l_state[g]:
            Tl[g] = prm["Th"]
        else:
            Tl[g] = prm["Tc"]
    t0 = 0
    print('\n\nProgress:')
    print('0%')
    tp = 0
    percent = 0
    for i in range(1, len(t_data)):
        if (i*prm["dt"])-t0 >= t_switch:
            l_state[n_switch%prm["Ngd"]] = not l_state[n_switch%prm["Ngd"]]
            if l_state[n_switch%prm["Ngd"]]:
                Tl[n_switch%prm["Ngd"]] = prm["Th"]
            else:
                Tl[n_switch%prm["Ngd"]] = prm["Tc"]
            n_switch += 1
            t0 = i*prm["dt"]
        if (i*prm["dt"])-tp >= prm["t"]/20:
            percent += 5
            tp = i*prm["dt"]
            print(f'{percent}%')
        Bt_data[i, :], BE_data[i, :], Bm_data[i, :] = gd_b_field(Bm_data[i-1, :], T_data[i-1, :], I_data[i-1, :])
        for g in range(prm["Ngd"]):
            T_data[i, g] = gd_temp_step(Bt_data[i-1, g], T_data[i-1, g], Bt_data[i, g], Tl[g])
        for g in range(prm["Ngd"]):
            U_data[i, g] = coil_volt(BE_data[i-1, g], BE_data[i, g])
            if i*prm["dt"]>1.0:
                I_data[i, g], q_data[i, g] = circuit_state(T_data[i, g], Bt_data[i, g], U_data[i, g], I_data[i-1, g], q_data[i-1, g])
                P_data[i, g] = gen_out(U_data[i, g], I_data[i, g])
    print('100%')
    results.append(np.array(main_values(t_data, T_data[:, 0], U_data[:, 0], I_data[:, 0], q_data[:, 0], P_data[:, 0], Bt_data[:, 0]), dtype=np.float64))
results = np.squeeze(np.array(results))

# single graphs

if not prm["series"] and prm["graphs"]:
    clrl = ['orangered', 'forestgreen', 'royalblue']
    
    for p in range(prm["Ngd"]):
        plt.plot(t_data[int(1/(prm["f"]*prm["dt"])):], T_data[int(1/(prm["f"]*prm["dt"])):, p], color=clrl[p], alpha=0.6)
    plt.axline((0, 295), (1000, 295), color="k", linestyle="--")
    plt.title("Gadolinium Temperature\n$T_{min} = " + f"{round(results[1], 2)}" + "$K, $T_{max} = " + f"{round(results[0], 2)}$K")
    plt.xlabel("Time $[s]$")
    plt.ylabel("Temperature $[K]$")
    plt.grid()
    plt.xlim(prm["t"]-1.0/prm["f"], prm["t"])
    plt.savefig("figs/tTs.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    for p in range(prm["Ngd"]):
        plt.plot(t_data[int(1/(prm["f"]*prm["dt"])):], Bt_data[int(1/(prm["f"]*prm["dt"])):, p], color=clrl[p], alpha=0.6)
    plt.title("Gadolinium Magnetic Induction")
    plt.xlabel("Time $[s]$")
    plt.ylabel("Magnetic induction $[T]$")
    plt.grid()
    plt.xlim(prm["t"]-1.0/prm["f"], prm["t"])
    plt.savefig("figs/tBs.pdf", format="pdf", bbox_inches="tight")
    plt.show() 
    
    for p in range(prm["Ngd"]):
        plt.plot(t_data[int(1/(prm["f"]*prm["dt"])):], U_data[int(1/(prm["f"]*prm["dt"])):, p], color=clrl[p], alpha=0.6)
    plt.title("Coil Electromotive Force")
    plt.xlabel("Time $[s]$")
    plt.ylabel("Electromotive force $[V]$")
    plt.grid()
    plt.xlim(prm["t"]-1.0/prm["f"], prm["t"])
    plt.savefig("figs/tVs.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    for p in range(prm["Ngd"]):
        plt.plot(t_data[int(1/(prm["f"]*prm["dt"])):], I_data[int(1/(prm["f"]*prm["dt"])):, p], color=clrl[p], alpha=0.6)
    plt.title("Output Current")
    plt.xlabel("Time $[s]$")
    plt.ylabel("Current $[A]$")
    plt.grid()
    plt.xlim(prm["t"]-1.0/prm["f"], prm["t"])
    plt.savefig("figs/tIs.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    for p in range(prm["Ngd"]):
        plt.plot(t_data[int(1/(prm["f"]*prm["dt"])):], q_data[int(1/(prm["f"]*prm["dt"])):, p], color=clrl[p], alpha=0.6)
    plt.title("Capacitor Electric Charge")
    plt.xlabel("Time $[s]$")
    plt.ylabel("Charge $[C]$")
    plt.grid()
    plt.xlim(prm["t"]-1.0/prm["f"], prm["t"])
    plt.savefig("figs/tQs.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    for p in range(prm["Ngd"]):
        plt.plot(t_data[int(1/(prm["f"]*prm["dt"])):], P_data[int(1/(prm["f"]*prm["dt"])):, p], color=clrl[p], alpha=0.6)
    plt.title("Output Power\n$U_{eff} = " + f"{round(1000*results[2], 1)}" + "$mV, $P_{eff} = " + f"{round(1000000000*results[3])}$nW")
    plt.xlabel("Time $[s]$")
    plt.ylabel("Power $[W]$")
    plt.grid()
    plt.xlim(prm["t"]-1.0/prm["f"], prm["t"])
    plt.savefig("figs/tPs.pdf", format="pdf", bbox_inches="tight")
    plt.show() 
    
    if prm["long_graphs"]:
        
        for p in range(prm["Ngd"]):
            plt.plot(t_data[int(1/prm["dt"]):], T_data[int(1/prm["dt"]):, p], color=clrl[p], alpha=0.6)
        plt.axline((0, 295), (1000, 295), color="k", linestyle="--")
        plt.title("Gadolinium Temperature")
        plt.xlabel("Time $[s]$")
        plt.ylabel("Temperature $[K]$")
        plt.grid()
        plt.xlim(1.0, prm["t"])
        plt.savefig("figs/tTl.pdf", format="pdf", bbox_inches="tight")
        plt.show()
        
        for p in range(prm["Ngd"]):
            plt.plot(t_data[int(1/prm["dt"]):], Bt_data[int(1/prm["dt"]):, p], color=clrl[p], alpha=0.6)
        plt.title("Gadolinium Magnetic Induction")
        plt.xlabel("Time $[s]$")
        plt.ylabel("Magnetic induction $[T]$")
        plt.grid()
        plt.xlim(1.0, prm["t"])
        plt.savefig("figs/tBl.pdf", format="pdf", bbox_inches="tight")
        plt.show() 
        
        for p in range(prm["Ngd"]):
            plt.plot(t_data[int(1/prm["dt"]):], U_data[int(1/prm["dt"]):, p], color=clrl[p], alpha=0.6)
        plt.title("Coil Electromotive Force")
        plt.xlabel("Time $[s]$")
        plt.ylabel("Electromotive force $[V]$")
        plt.grid()
        plt.xlim(1.0, prm["t"])
        plt.savefig("figs/tVl.pdf", format="pdf", bbox_inches="tight")
        plt.show()
        
        for p in range(prm["Ngd"]):
            plt.plot(t_data[int(1/prm["dt"]):], I_data[int(1/prm["dt"]):, p], color=clrl[p], alpha=0.6)
        plt.title("Output Current")
        plt.xlabel("Time $[s]$")
        plt.ylabel("Current $[A]$")
        plt.grid()
        plt.xlim(1.0, prm["t"])
        plt.savefig("figs/tIl.pdf", format="pdf", bbox_inches="tight")
        plt.show()
        
        for p in range(prm["Ngd"]):
            plt.plot(t_data[int(1/prm["dt"]):], q_data[int(1/prm["dt"]):, p], color=clrl[p], alpha=0.6)
        plt.title("Capacitor Electric Charge")
        plt.xlabel("Time $[s]$")
        plt.ylabel("Charge $[C]$")
        plt.grid()
        plt.xlim(1.0, prm["t"])
        plt.savefig("figs/tQl.pdf", format="pdf", bbox_inches="tight")
        plt.show()
        
        for p in range(prm["Ngd"]):
            plt.plot(t_data[int(1/prm["dt"]):], P_data[int(1/prm["dt"]):, p], color=clrl[p], alpha=0.6)
        plt.title("Output Power")
        plt.xlabel("Time $[s]$")
        plt.ylabel("Power $[W]$")
        plt.grid()
        plt.xlim(1.0, prm["t"])
        plt.savefig("figs/tPl.pdf", format="pdf", bbox_inches="tight")
        plt.show() 

# series graphs -- manual change of x-axis is required

if prm["series"] and prm["graphs"]:
    
    plt.plot(range_x, results[:, 0], color = 'darkred')
    plt.plot(range_x, results[:, 1], color = 'darkblue')
    plt.fill_between(range_x, results[:, 0], results[:, 1], color='lavender')
    plt.axline((range_x[0], 295), (range_x[-1], 295), color="k", linestyle="--")
    plt.title("Temperature Range")
    plt.xlabel(series_x)
    plt.ylabel("Temperature $[K]$")
    plt.grid()
    plt.savefig("figs/xTr.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    plt.plot(range_x, results[:, 2], color = 'deepskyblue')
    plt.title("Effective Voltage")
    plt.xlabel(series_x)
    plt.ylabel("Voltage $[V]$")
    plt.grid()
    plt.savefig("figs/xVr.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    plt.plot(range_x, results[:, 3], color = 'darkorange')
    plt.title("Effective Power")
    plt.xlabel(series_x)
    plt.ylabel("Power $[W]$")
    plt.grid()
    plt.savefig("figs/xPr.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    plt.plot(range_x, results[:, 4], color = 'deeppink')
    plt.title("Average Gadolinium Magnetic Induction")
    plt.xlabel(series_x)
    plt.ylabel("Magnetic induction $[T]$")
    plt.grid()
    plt.savefig("figs/xBr.pdf", format="pdf", bbox_inches="tight")
    plt.show()

