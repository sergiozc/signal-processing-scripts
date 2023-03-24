# -*- coding: utf-8 -*-
"""

@author: Sergio

"""

import numpy as np
import scipy
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.io import wavfile
from scipy import signal
import wave
import sounddevice as sd


def calculo_snr(noise, power):
    return 10*np.log10(power/noise)


def directividad_teo(barrido_exp, barrido_polar, f, n, d, c, w):
    
    '''Función que calcula la directividad teórica para una frecuencia concreta y 
    se proporciona el valor de la directividad polar para ser representada en 
    coordenadas polares'''
    
    D = np.zeros(len(barrido_exp), dtype = 'complex_')
    Dtot = np.zeros(len(barrido_polar), dtype = 'complex')
    
    #Fórmula página 15 transparencias
    for i in range (len(barrido_exp)):
        Dsum = 0
        for n in range (N):
            Dsum = Dsum + np.exp(2 * np.pi * 1j * f * n * d * np.cos(barrido_exp[i]) / c) * np.conj(w[n])
    
        D[i] = Dsum

    #Parte simétrica de la directividad
    Dsim = D[::-1]
    Dtot = np.concatenate((D, Dsim), axis=0)
    
    return np.abs(Dtot)

def calcula_pesos(N, d, phi, c, f):
    
    '''Función que calcula los pesos del beamformer delay and sum dado un ángulo
    de incidencia (phi), la distancia, la velocidad de propagación y la frecuencia'''
    
    w = np.zeros((N), dtype = 'complex_')   # Pesos beamformer (complejos)
    n = np.arange(N) 
    retardo = (n * d * np.cos(phi)) / c     # Caso array linear uniforme

    #Pesos
    for i in range (N):
        w[i] = np.exp(1j * 2 * np.pi * f * retardo[i])
    
    w = w / N   # Contribución uniforme
    
    return w
    

def directividad_frecuencia(f_barrido, barrido_exp, N, d, phi, c):
    
    ''' Función que calcula la directividad del beamformer en función del ángulo
    theta y en función de la frecuencia'''
    
    Df = np.zeros([len(f_barrido), len(barrido_exp)], dtype = 'complex_')
    
    for f in range(len(f_barrido)):
        # Calculamos los pesos correspondientes a todas las frecuencias
        wf = calcula_pesos(N, d, phi, c, f_barrido[f])
        for theta in range (len(barrido_exp)):
            Dsum = 0
            for n in range (N):
                Dsum = Dsum + np.exp(2 * np.pi * 1j * f_barrido[f] * n * d * np.cos(barrido_exp[theta]) / c) * np.conj(wf[n])
    
            Df[f, theta] = Dsum
    
    return np.abs(Df)


#%% Parámetros del Beamformer

plt.close('all')

load = sio.loadmat('signals_array.mat')
objects_signals = load['xc']
signals = np.zeros([44057, 7])
signals = np.column_stack((objects_signals[0, 0], objects_signals[0 ,1], objects_signals[0, 2], \
objects_signals[0 ,3], objects_signals[0, 4], objects_signals[0 ,5], objects_signals[0, 6]))

# wavfile.write("signal_orig.wav", 44100, signals[:,1])


phi = np.pi / 2     # Ángulo de incidencia de la fuente
f = 16000           # Frecuencia de muestreo
c = 340             # Velocidad del sonido
N = 7               # Número de elementos del array
d = 0.04            # Separación de los elementos del array
n = np.arange(N)    # Canales
w = np.zeros((N), dtype = 'complex_')   # Pesos beamformer (complejos)
y = np.zeros(len(signals))  # Señal a la salida del beamformer

barrido_exp = np.linspace(0, np.pi, 200)     #Barrido experimental (0-180 grados)
barrido_polar = np.linspace(0, 2*np.pi, 400) #Barrido para representar (360 grados)
f_barrido = np.linspace(100, 8000, 200)     # Barrido en frecuencias de 100 a 8000 Hz
wf = np.zeros(N, dtype = 'complex_')
Df = np.zeros([len(f_barrido), len(barrido_exp)], dtype = 'complex_')


#%% CÁLCULO DE LOS PESOS DEL BEAMFORMER

retardo = (n * d * np.cos(phi)) / c     # Caso array linear uniforme

#Pesos
for i in range (N):
    w[i] = np.exp(1j * 2 * np.pi * f * retardo[i])

w = w / N     # Contribución uniforme

#%% PROMEDIO (sum)

for j in range (N):
    y = y + (signals[:, j] * w[j])
    
plt.figure(14)
plt.plot(y)
plt.title('Señal tras "SUM"')
plt.xlabel('Muestras')
plt.ylabel('Amplitud')

# wavfile.write("signal.wav", f,np.abs(y))

#%% CÁLCULO DE LA SNR

noise_before = np.var(signals[0:3000, 3])
power_before = np.var(signals[3001:, 3])
snr_before = calculo_snr(noise_before, power_before)
print("SNR antes: " + str(snr_before))

noise_after = np.var(y[0:3000])
power_after = np.var(y[3001:])
snr_after = calculo_snr(noise_after, power_after)
print("SNR después: " + str(snr_after))


#%% PATRÓN DIRECTIVIDAD TEÓRICO

Dteo = directividad_teo(barrido_exp, barrido_polar, f, n, d, c, w)

plt.figure(15)
plt.polar(barrido_polar, abs(Dteo))
plt.title('Diagrama directividad teórica')


#%% PATRÓN DE DIRECTIVIDAD EN FRECUENCIAS

Df = directividad_frecuencia(f_barrido, barrido_exp, N, d, phi, c)


plt.figure(17)
plt.pcolormesh(f_barrido, barrido_exp, np.abs(Df.T), shading = 'auto')
plt.colorbar()
plt.xlabel('Frecuencia (Hz)')
plt.ylabel('θ (rad)')
plt.title('Directividad en función de la frecuencia y ángulo')
plt.show()



# Representación de la directividad a una frecuencia concreta para
# distintas distancias entre elementos del array
dist = [0.02, 0.04, 0.1]
plt.figure(16)
for dis in dist:
    Df = directividad_frecuencia(f_barrido, barrido_exp, N, dis, phi, c)
    plt.plot(barrido_exp, np.abs(Df[50, :]), label= 'd = ' + str(dis))

plt.legend()
plt.xlabel('θ (rad)')
plt.ylabel('|D(f, θ)|')
plt.title('Directividad en función de theta para distintas distancias')



