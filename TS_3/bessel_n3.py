#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 17:43:40 2026

@author: lrobledo
"""

# Librerías externas NumPy, SciPy y Matplotlib
from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np 
from scipy.signal import TransferFunction, bode
from scipy.signal import group_delay

# Librería de TC2, esta la vas a usar mucho
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot

from sympy import *

#transferencia bessel orden 3
num = [15]   # 1
den = [1, 6, 15, 15]
#raices
np.roots(den)

#verificacion salley key
# Variables simbólicas
V1, Vo, Vx, s, C1, R1, R2, C2 = symbols('V1 Vo Vx s C1 R1 R2 C2')
# Ecuaciones
eq1 = Eq((V1-Vx)/R1, (Vx-Vo)/R2+(Vx-Vo)*s*C1)
eq2 = Eq((Vx-Vo)/R2, Vo*s*C2)
# Resolver sistema
sol = solve((eq1, eq2), (Vo, Vx))[Vo]/V1
print(sol)

Dd = ( 2.3222 ) / ( 2.3222**2 + 2.5**2 ) + ( 1.8389 ) / ( 1.8389**2  + ( 2.5**2 - 1.7544 )**2 ) + ( 1.8389 ) / ( 1.8389**2  + ( 2.5**2 + 1.7544 )**2 )

print_latex(Dd)

my_tf = TransferFunction(num, den)
sys = TransferFunction(num, den)
plt.close('all')

plt.show()
bodePlot(my_tf, fig_id=1, filter_description = 'Cn=5')

pzmap(my_tf, fig_id=2, filter_description = 'Cn=5')


# cálculo del retardo de grupo
w, gd = group_delay((num, den), w=10000)

# frecuencia de interés
w_obj = 2.5

# índice más cercano
idx = np.argmin(np.abs(w - w_obj))

# valor
gd_obj = gd[idx]

print(f"Retardo de grupo en {w[idx]:.3f} rad/s = {gd_obj:.5f}")



# w = np.logspace(0, 0.43, 100000)

# w, mag, phase = bode(sys, w=w)

# fig, ax = plt.subplots(2, 1)
# # Magnitud
# ax[0].semilogx(w, mag)
# ax[0].grid(True)
# ax[0].set_ylabel('Magnitud [dB]')

# # Fase
# ax[1].semilogx(w, phase*np.pi/180)
# ax[1].grid(True)
# ax[1].set_ylabel('Fase [rad]')
# ax[1].set_xlabel('Frecuencia angular [rad/s]')