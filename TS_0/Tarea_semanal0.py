#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 18:35:37 2026

@author: lrobledo
"""

# Librerías externas NumPy, SciPy y Matplotlib
from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np


# Librería de TC2, esta la vas a usar mucho
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot

#componentes
r1 = 10
r2= 10
r3= 10
c1=1/10
k = 1/(c1*r3)

num = [1, -k*r2/r1]   # s - k*r1
den = [1, k]       # s + k

#my_tf = TransferFunction(num, den)

plt.close('all')
relaciones = [1, 0.5, 2]

for rel in relaciones:
    
    num = [1, -k*rel]   # s - k*(R2/R1)
    den = [1, k]        # s + k
    
    my_tf = TransferFunction(num, den)
    
    bodePlot(my_tf, fig_id=1, filter_description = 'R2/R1={:3.2f}'.format(rel))
    
    pzmap(my_tf, fig_id=2, filter_description = 'R2/R1={:3.2f}'.format(rel))
    
    GroupDelay(my_tf, fig_id=3, filter_description = 'R2/R1={:3.2f}'.format(rel))



