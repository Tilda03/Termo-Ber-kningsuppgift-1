import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, quad


C = 10
R = 296.8 #J/kg * K
P_c = 3.39 * 10**6 #Pascal
T_c = 126.2 #Kelvin

P_1 = 500 * 10**3 #Pascal

T_20 = 293 #Kelvin
T_140 = 133 #Kelvin

a = (27 * R**2 * T_c**2)/(64 * P_c)
b = (R * T_c)/(8 * P_c)


def compression_to_v(compression_factor, T):
    v1 = R * T / P_1
    v2 = v1 / compression_factor    
    return v1, v2

# Ideal gas

def pressure_ideal(v, T):
    return R * T / v


# Diderik vadn der Waals

def pressure_waals(v, T):
    return (R * T)/(v-b) - (a/v**2)


def work(P, compression_factor, P_1, T):
    v1 = R * T / P_1
    v2 = v1 / compression_factor
    I, err = quad(P, v1, v2, args=T)
    return I


v = np.linspace(compression_to_v(C, 293), 800)

print(work(pressure_ideal,C,P_1, 293))

print(work(pressure_waals, C, P_1, 293))

print(work(pressure_ideal,C,P_1, 133))


print(work(pressure_waals, C, P_1, 133))


plt.plot(v, pressure_ideal(v, 293), label='ideal gas')
plt.legend()
plt.show()

plt.plot(v, pressure_waals(v, 293), label='waals gas')

plt.legend()
plt.show()
