import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.optimize import brentq


T = int(input('Enter temperature (c): '))+273
C = int(input('Enter compression factor (%): '))
P_1 = int(input('Enter starting pressure (kPa): '))*10**3


R = 296.8 #J/kg * K
P_c = 3.39 * 10**6 #Pascal
T_c = 126.2 #Kelvin


a = (27 * R**2 * T_c**2)/(64 * P_c)
b = (R * T_c)/(8 * P_c)


def compression_to_v_ideal(P_1, R, compression_factor, T):
    v1 = R * T / P_1
    v2 = v1 / compression_factor 
    return np.linspace(v1, v2, 201)

def compression_to_v_waals(P_1, R, compression_factor, T):
    def f(v):
        return (P_1+a/v**2)*(v-b)-R*T
    
    #Vi vet att xl måste vara större än b
    xl = b * 1.001
    xr = 1.0

    v1 = brentq(f, xl, xr)
    v2 = v1 / compression_factor 

    return np.linspace(v1, v2, 201)

# Ideal gas

def pressure_ideal(R, v, T):
    return R * T / v


# Diderik vadn der Waals

def pressure_waals(R, v, T):
    return (R * T)/(v-b) - (a/v**2)


def work(P, v_range):
    I = trapezoid(P, v_range)
    return I


data = np.loadtxt("mättnadskurva_kväve.txt")
v_sat = data[:,2]
P_sat = data[:,0]

# Ranges for v and for p
v_ideal_range = compression_to_v_ideal(P_1, R, C, T)
P_ideal_range = pressure_ideal(R, v_ideal_range, T)

v_waals_range = compression_to_v_waals(P_1, R, C, T)
P_waals_range = pressure_waals(R, v_waals_range, T)

print(f'Total work for ideal gas law: {work(P_ideal_range, v_ideal_range):.0f} J')
print(f'Total work for waals gas law: {work(P_waals_range, v_waals_range):.0f} J')

# plot for ideal and waals
plt.subplot(1, 2, 1)
plt.plot(v_ideal_range, P_ideal_range, label='Ideal', color='royalblue')
plt.plot(v_waals_range, P_waals_range, label='Waals', color='crimson')
plt.xlabel('Specific volume [m^3/kg]')
plt.ylabel('Pressure [Pa]')
plt.legend()
plt.grid(True, which="both")

# plot for saturation curve
plt.subplot(1, 2, 2)
plt.plot(v_sat, P_sat, label='Saturation', color='forestgreen') 
plt.xlabel('Specific volume [m^3/s]')
plt.ylabel('Pressure [kPa]')
plt.legend()
plt.xscale('log')
plt.yscale('log')

plt.show()
