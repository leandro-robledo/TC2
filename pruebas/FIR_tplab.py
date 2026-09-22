#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 14:39:14 2026

@author: lrobledo
"""

# clase profe

import sympy as sp
import numpy as np
import scipy.signal as sig
from scipy.signal.windows import hamming, kaiser, blackmanharris
import matplotlib.pyplot as plt

from pytc2.sistemas_lineales import plot_plantilla, group_delay


# frecuencia de muestreo normalizada
fs = 2.0
# tamaño de la respuesta al impulso
cant_coef = 51
filter_type = 'lowpass'

fpass = 0.2 # 
ripple = 1 # dB
fstop = 0.6 # Hz
attenuation = 60 # dB
fneg= 0.8
segundpaso=20

# construyo la plantilla de requerimientos
frecs = [0.0,  fpass,     fstop, 1.0]
gains = [0,   -ripple, -attenuation, -np.inf] # dB

gains = 10**(np.array(gains)/20)

# algunas ventanas para evaluar
#win_name = 'boxcar'
#win_name = 
win_name = kaiser
#win_name = 'flattop'

# FIR design
num_bh = sig.firwin2(cant_coef, frecs, gains , window='blackmanharris' )
num_hm = sig.firwin2(cant_coef, frecs, gains , window='hamming' )
num_ka = sig.firwin2(cant_coef, frecs, gains , window=('kaiser',14))
den = 1.0

#num_firls= sig.firls(cant_coef,frecs,[1,1,0,0],fs=fs)

def plot_freq_resp_fir(this_num, this_desc):

    wrad, hh = sig.freqz(this_num, 1.0)
    ww = wrad / np.pi
    
    plt.figure(1)

    plt.plot(ww, 20 * np.log10(abs(hh)), label=this_desc)

    plt.title('FIR diseñado por métodos directos - Taps:' + str(cant_coef) )
    plt.xlabel('Frequencia normalizada')
    plt.ylabel('Modulo [dB]')
    plt.grid(which='both', axis='both')

    axes_hdl = plt.gca()
    axes_hdl.legend()
    
    plt.figure(2)

    phase = np.unwrap(np.angle(hh))

    plt.plot(ww, phase, label=this_desc)

    plt.title('FIR diseñado por métodos directos - Taps:' + str(cant_coef))
    plt.xlabel('Frequencia normalizada')
    plt.ylabel('Fase [rad]')
    plt.grid(which='both', axis='both')

    axes_hdl = plt.gca()
    axes_hdl.legend()

    plt.figure(3)

    # ojo al escalar Omega y luego calcular la derivada.
    gd_win = group_delay(wrad, phase)

    plt.plot(ww, gd_win, label=this_desc)

    plt.ylim((np.min(gd_win[2:-2])-1, np.max(gd_win[2:-2])+1))
    plt.title('FIR diseñado por métodos directos - Taps:' + str(cant_coef))
    plt.xlabel('Frequencia normalizada')
    plt.ylabel('Retardo [# muestras]')
    plt.grid(which='both', axis='both')

    axes_hdl = plt.gca()
    axes_hdl.legend()    
plt.close("all")
plot_freq_resp_fir(num_bh, filter_type+ '-blackmanharris')    
plot_freq_resp_fir(num_hm, filter_type+ '-hamming')    
plot_freq_resp_fir(num_ka, filter_type+ '-kaiser-b14')    

#plot_freq_resp_fir(num_firls, filter_type+ 'firls')    
    
# sobreimprimimos la plantilla del filtro requerido para mejorar la visualización    
fig = plt.figure(1)    
plot_plantilla(filter_type = filter_type , fpass = fpass, ripple = ripple , fstop = fstop, attenuation = attenuation, fs = fs)
ax = plt.gca()
ax.legend()

# reordenamos las figuras en el orden habitual: módulo-fase-retardo
plt.figure(2)    
axes_hdl = plt.gca()
axes_hdl.legend()

plt.figure(3)    
axes_hdl = plt.gca()
axes_hdl.legend()

plt.show()



#%%
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

from pytc2.sistemas_lineales import plot_plantilla, group_delay


# ============================================================
# ESPECIFICACIONES
# ============================================================

fs = 2.0                  # frecuencia de muestreo normalizada
cant_coef = 21           # cantidad de coeficientes FIR

filter_type = 'lowpass'

fpass = 0.2               # fin de banda pasante
fstop = 0.6               # comienzo de banda de rechazo

ripple = 1                # ripple permitido en dB
attenuation = 60          # atenuación mínima en dB


# ============================================================
# PLANTILLA
# ============================================================

frecs = [0.0, fpass, fstop, fs/2]

gains = [0, -ripple, -attenuation, -attenuation]


# ============================================================
# FIR POR REMEZ / PARKS-MCCLELLAN
# ============================================================

num_remez = sig.remez(
    cant_coef,
    [0, fpass, fstop, fs/2],
    [1, 0],
    weight=[1, 1],
    fs=fs
)


# ============================================================
# FUNCIÓN PARA GRAFICAR
# ============================================================

def plot_freq_resp_fir(this_num, this_desc):

    wrad, hh = sig.freqz(this_num, 1.0)

    # frecuencia normalizada respecto de Nyquist
    ww = wrad / np.pi

    # -------------------------
    # MÓDULO
    # -------------------------

    plt.figure(1)

    modulo_db = 20 * np.log10(np.maximum(np.abs(hh), 1e-12))

    plt.plot(
        ww,
        modulo_db,
        label=this_desc
    )

    plt.title(
        'FIR - Taps: ' + str(cant_coef)
    )

    plt.xlabel('Frecuencia normalizada')
    plt.ylabel('Módulo [dB]')

    plt.grid(which='both', axis='both')

    plt.legend()


    # -------------------------
    # FASE
    # -------------------------

    plt.figure(2)

    phase = np.unwrap(np.angle(hh))

    plt.plot(
        ww,
        phase,
        label=this_desc
    )

    plt.title(
        'FIR - Taps: ' + str(cant_coef)
    )

    plt.xlabel('Frecuencia normalizada')
    plt.ylabel('Fase [rad]')

    plt.grid(which='both', axis='both')

    plt.legend()


    # -------------------------
    # RETARDO DE GRUPO
    # -------------------------

    plt.figure(3)

    gd = group_delay(wrad, phase)

    plt.plot(
        ww,
        gd,
        label=this_desc
    )

    plt.title(
        'FIR - Taps: ' + str(cant_coef)
    )

    plt.xlabel('Frecuencia normalizada')
    plt.ylabel('Retardo [muestras]')

    plt.grid(which='both', axis='both')

    plt.legend()


# ============================================================
# GRAFICAR FILTRO
# ============================================================

plt.close("all")

plot_freq_resp_fir(
    num_remez,
    filter_type + '-Remez'
)


# ============================================================
# SUPERPONER PLANTILLA
# ============================================================

plt.figure(1)

plot_plantilla(
    filter_type=filter_type,
    fpass=fpass,
    ripple=ripple,
    fstop=fstop,
    attenuation=attenuation,
    fs=fs
)

plt.legend()


plt.figure(2)
plt.legend()

plt.figure(3)
plt.legend()

plt.show()


#%% prueba mia xd
from pytc2.filtros_digitales import fir_design_ls

import sympy as sp
import numpy as np
import scipy.signal as sig
from scipy.signal.windows import hamming, kaiser, blackmanharris
import matplotlib.pyplot as plt

from pytc2.sistemas_lineales import plot_plantilla, group_delay


Ftype = 'm'
# frecuencia de muestreo normalizada
fs = 2.0
# tamaño de la respuesta al impulso
cant_coef = 51
filter_type = 'lowpass'

fpass = 0.2 # 
ripple = 1 # dB
fstop = 0.5 # Hz
attenuation = 60 # dB
fneg= 0.8
segundpaso=20

# construyo la plantilla de requerimientos
frecs = [0.0,  fpass,     fstop,          1.0]
gains = [0,   -ripple, -attenuation,  0] # dB

gains = 10**(np.array(gains)/20)

Be = [0, 0.4, 0.41, 1.0]
D = [1, 1, 0, 0]

Be_sig = Be
D_sig = [1., 0.]

# enfatizamos stop
W = [1., 50.]

num_firls= sig.firls(cant_coef,frecs,[1,1,0,0],fs=fs)

#hh_mi_ls_stop = fir_design_ls(cant_coef ,band_edges=Be, desired=D, 
 #                         weight=W, filter_type = Ftype, 
  #                        grid_density= 4)

# enfatizamos ripple en paso
W = [50., 1.]

#hh_mi_ls_paso = fir_design_ls(order=N, band_edges=Be, desired=D, 
#                          weight=W, filter_type = Ftype, 
#                          grid_density= 4)

#H_mi_ls_stop = np.fft.fft(hh_mi_ls_stop, fft_sz)
#H_mi_ls_paso = np.fft.fft(hh_mi_ls_paso, fft_sz)

def plot_freq_resp_fir(this_num, this_desc):

    wrad, hh = sig.freqz(this_num, 1.0)
    ww = wrad / np.pi
    
    plt.figure(1)

    plt.plot(ww, 20 * np.log10(abs(hh)), label=this_desc)

    plt.title('FIR diseñado por métodos directos - Taps:' + str(cant_coef) )
    plt.xlabel('Frequencia normalizada')
    plt.ylabel('Modulo [dB]')
    plt.grid(which='both', axis='both')

    axes_hdl = plt.gca()
    axes_hdl.legend()
    
    plt.figure(2)

    phase = np.unwrap(np.angle(hh))

    plt.plot(ww, phase, label=this_desc)

    plt.title('FIR diseñado por métodos directos - Taps:' + str(cant_coef))
    plt.xlabel('Frequencia normalizada')
    plt.ylabel('Fase [rad]')
    plt.grid(which='both', axis='both')

    axes_hdl = plt.gca()
    axes_hdl.legend()

    plt.figure(3)

    # ojo al escalar Omega y luego calcular la derivada.
    gd_win = group_delay(wrad, phase)

    plt.plot(ww, gd_win, label=this_desc)

    plt.ylim((np.min(gd_win[2:-2])-1, np.max(gd_win[2:-2])+1))
    plt.title('FIR diseñado por métodos directos - Taps:' + str(cant_coef))
    plt.xlabel('Frequencia normalizada')
    plt.ylabel('Retardo [# muestras]')
    plt.grid(which='both', axis='both')

    axes_hdl = plt.gca()
    axes_hdl.legend() 



# Graficar la respuesta en frecuencia
#plt.plot(frecs, 20*np.log10(np.abs(H[:hh_mi_ls_stop])), label='LS')
#plt.plot(frecuencias[:half_fft_sz], 20*np.log10(np.abs(H_mi_ls_paso[:half_fft_sz])), label='paso')
#plt.plot(frecuencias[:half_fft_sz], 20*np.log10(np.abs(H_mi_ls_stop[:half_fft_sz])), label='stop')
num_firls= sig.firls(cant_coef,frecs,[1,1,0,0],fs=fs)
plot_freq_resp_fir(num_firls, filter_type+ 'firls') 
 
plt.title("Respuesta en Frecuencia del Filtro FIR Diseñado")
plt.xlabel("Frecuencia Normalizada")
plt.ylabel("Magnitud")
plt.legend()

# sobreimprimimos la plantilla del filtro requerido para mejorar la visualización    
fig = plt.figure(1)    
plot_plantilla(filter_type = filter_type , fpass = fpass, ripple = ripple , fstop = fstop, attenuation = attenuation, fs = fs)
ax = plt.gca()
ax.legend()

# reordenamos las figuras en el orden habitual: módulo-fase-retardo
plt.figure(2)    
axes_hdl = plt.gca()
axes_hdl.legend()

plt.figure(3)    
axes_hdl = plt.gca()
axes_hdl.legend()

plt.show()


