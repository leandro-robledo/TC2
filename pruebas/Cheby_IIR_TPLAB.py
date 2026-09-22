#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 22:32:20 2026

@author: lrobledo
"""

import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt


# ============================================================
# ESPECIFICACIONES DEL FILTRO
# ============================================================

fs = 20000.0          # Frecuencia de muestreo [Hz]
f_nyquist = fs / 2

fp = 100.0           # Fin de banda pasante [Hz]
fstop = 300.0        # Inicio de banda de rechazo [Hz]

gpass = 1.0          # Máxima atenuación banda pasante [dB]
gstop = 60.0         # Mínima atenuación banda de rechazo [dB]


# ============================================================
# DISEÑO DEL FILTRO CHEBYSHEV TIPO I
# ============================================================

sos = sig.iirdesign(
    wp=fp,
    ws=fstop,
    gpass=gpass,
    gstop=gstop,
    analog=False,
    ftype='cheby1',
    output='sos',
    fs=fs
)


# ============================================================
# INFORMACIÓN DEL FILTRO
# ============================================================

cantidad_sos = sos.shape[0]
orden = 2 * cantidad_sos

print("==============================================")
print(" FILTRO CHEBYSHEV TIPO I")
print("==============================================")

print(f"Frecuencia de muestreo = {fs:.1f} Hz")
print(f"Frecuencia de Nyquist   = {f_nyquist:.1f} Hz")
print(f"Frecuencia pasante      = {fp:.1f} Hz")
print(f"Frecuencia de rechazo   = {fstop:.1f} Hz")
print(f"Ripple permitido        = {gpass:.1f} dB")
print(f"Atenuación mínima       = {gstop:.1f} dB")
print(f"Orden del filtro        = {orden}")
print(f"Cantidad de SOS         = {cantidad_sos}")

print("\nCoeficientes SOS:")
print(sos)


# ============================================================
# RESPUESTA EN FRECUENCIA
# ============================================================

w, h = sig.sosfreqz(
    sos,
    worN=4096,
    fs=fs
)

H_dB = 20 * np.log10(np.maximum(np.abs(h), 1e-12))


# ============================================================
# FRECUENCIAS NORMALIZADAS
# ============================================================

w_normalizada = w / f_nyquist

fp_normalizada = fp / f_nyquist
fstop_normalizada = fstop / f_nyquist


# ============================================================
# GRÁFICA
# ============================================================

plt.figure(figsize=(10, 6))

plt.plot(
    w_normalizada,
    H_dB,
    label='Chebyshev tipo I'
)

# Límites de la banda pasante
plt.axvline(
    fp_normalizada,
    linestyle='--',
    label='fp = 100 Hz'
)

# Inicio banda de rechazo
plt.axvline(
    fstop_normalizada,
    linestyle='--',
    label='fstop = 300 Hz'
)

# Límite de ripple
plt.axhline(
    -gpass,
    linestyle=':',
    label='-1 dB'
)

# Límite de atenuación
plt.axhline(
    -gstop,
    linestyle=':',
    label='-60 dB'
)


plt.title('Filtro IIR - Chebyshev tipo I')
plt.xlabel('Frecuencia normalizada respecto de Nyquist')
plt.ylabel('Magnitud [dB]')

plt.xlim(0, 1)
plt.ylim(-100, 5)

plt.grid(True)
plt.legend()

plt.show()