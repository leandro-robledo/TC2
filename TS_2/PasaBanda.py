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


# Librería de TC2, esta la vas a usar mucho
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot

from sympy import *

# Variables simbólicas
V1, Vo, Vx, s, C, R1, R2, R3 = symbols('V1 Vo Vx s C R1 R2 R3')

# Ecuaciones
eq1 = Eq(V1/R1, -Vo*(s*C + 1/R2)-Vx/R3)
eq2 = Eq(Vo/R3, Vx*s*C)

# Resolver sistema
sol = solve((eq1, eq2), (Vo, Vx))[Vo]/V1

print(sol)

