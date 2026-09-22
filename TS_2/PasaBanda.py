#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 19:34:17 2026

@author: lrobledo
"""

# Librerías externas NumPy, SciPy y Matplotlib
from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np
from pytc2.sistemas_lineales import pretty_print_lti

from pytc2.sistemas_lineales import analyze_sys, pretty_print_bicuad_omegayq, tf2sos_analog, pretty_print_SOS

from pytc2.general import print_subtitle
# Librería de TC2, esta la vas a usar mucho
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot

from sympy import *

# Variables simbólicas
V1, Vo, Vx, Vy, s, C, R, R1, R3, R4, R5  = symbols('V1 Vo Vx Vy s C R R1 R3 R4 R5')

# Ecuaciones
eq1 = Eq((V1-Vx)/R, Vx*s*C+(Vx-Vy)/R1)
eq2 = Eq((Vy-Vx)*-s*C, (Vx-Vo)/R3)
eq3 = Eq((Vo-Vx)*R4, Vx/R5)


# Resolver sistema
sol = solve((eq1, eq2, eq3), (Vo, Vx, Vy), dict=True)

print (sol)

#no es el mismo
#parcial del curso de los lunes

Q=7.975
psy= 0.122
W=140016.9
#desnormalizada
#num = [(W**2)/(Q**2),0,0]   # 1
#den = [1,1.41*W/(Q),(2+1/Q**2)*W**2, (1.41*W**3)/(Q),W**4]       # s + k

#normalizada
num = [(1)/(Q**2),0,0]   # 1
den = [1,1.41/(Q),(2+1/Q**2), (1.41)/(Q),1]


my_tf = TransferFunction(num, den)

plt.close('all')

my_tf = TransferFunction(num, den)

bodePlot(my_tf, fig_id=1, filter_description = 'Cn=5')

pzmap(my_tf, fig_id=2, filter_description = 'Cn=5')

GroupDelay(my_tf, fig_id=3, filter_description = 'Cn=5')

this_sos = tf2sos_analog(num, den)
pretty_print_SOS(this_sos, mode='omegayq')
pretty_print_lti(num,den)