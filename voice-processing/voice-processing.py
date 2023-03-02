# -*- coding: utf-8 -*-
"""
Modelado y procesamiento de voz

Autores: Sergio Zapata Caparrós y Antonio Simón Martín

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from spectrum import *
from math import dist


def vocal_x(indice):
    """
    Función que a partir de un segmento vocálico
    devuelve el cepstrum_LPC correspondiente
    Argumento de entrada: indice segmento
    Return: cepstrum_LPC
    """
    if indice == 1:
        x = np.loadtxt('vocales/voc_a.asc')
    elif indice == 2:
        x = np.loadtxt('vocales/voc_e.asc')
    elif indice == 3:
        x = np.loadtxt('vocales/voc_i.asc')
    elif indice == 4:
        x = np.loadtxt('vocales/voc_o.asc')
    elif indice == 5:
        x = np.loadtxt('vocales/voc_u.asc')
    else:
        x = np.loadtxt('vocales/voc_x.asc')
        
    p = 12 #número de polos para el filtro
    a = np.empty(12)

    a_,sig2,k = aryule(x,p)
    a = np.concatenate(([1],a_), axis = 0)
    #Calculamos el espectro AR con la fórmula teórica
    Px_AR = sig2 / (np.abs(np.fft.fft(a,Lfft)) ** 2)
    
    cepstrum_LPC = np.fft.ifft(np.log(Px_AR))
    
    cepstrum = cepstrum_LPC[1:12]
    
    c = np.empty(len(cepstrum))
    
    for i in range (len(cepstrum)):
        c[i]=abs(cepstrum[i])
    
    return c

def distancia_minima (cepstrum_a, cepstrum_e, cepstrum_i, cepstrum_o, cepstrum_u, cepstrum_x):
    """
    Función que calcula la distancia mínima al cepstrum_x para hallar la vocal incógnita
    El cepstrum que esté más cerca, corresponderá al cepstrum de la vocal oculta.
    Argumento de entrada: cepstrum de cada vocal + la incógnita
    Return: indice de la vocal correspondiente al cepstrum más cercano
    """
    distancia_a = np.empty(len(cepstrum_x))
    distancia_e = np.empty(len(cepstrum_x))
    distancia_i = np.empty(len(cepstrum_x))
    distancia_o = np.empty(len(cepstrum_x))
    distancia_u = np.empty(len(cepstrum_x))
    
    for i in range (len(cepstrum_x)):
        distancia_a[i] = (np.abs(cepstrum_x[i] - cepstrum_a[i])) ** 2
        distancia_e[i] = (np.abs(cepstrum_x[i] - cepstrum_e[i])) ** 2
        distancia_i[i] = (np.abs(cepstrum_x[i] - cepstrum_i[i])) ** 2
        distancia_o[i] = (np.abs(cepstrum_x[i] - cepstrum_o[i])) ** 2
        distancia_u[i] = (np.abs(cepstrum_x[i] - cepstrum_u[i])) ** 2
        
    media_a = np.mean(distancia_a)
    media_e = np.mean(distancia_e)
    media_i = np.mean(distancia_i)
    media_o = np.mean(distancia_o)
    media_u = np.mean(distancia_u)
        
    distancia = np.array([media_a, media_e, media_i, media_o, media_u])
        
    indice_vocal = np.argmin(distancia)
        
    
    return indice_vocal


def identifica_vocal(indice):
    """
    Dado el índice de la mínima distancia del cepstrum,
    identifica la vocal correspondiente
    """
    if indice == 0:
        print('VOCAL A')
    elif indice == 1:
        print('VOCAL E')
    elif indice == 2:
        print('VOCAL I')
    elif indice == 3:
        print('VOCAL O')
    else:
        print('VOCAL U')
    pass
    
#%% ANÁLISIS VOCAL "E"

Lfft = 1024
N = 256

e = np.loadtxt('vocales/voc_e.asc')
plt.figure(1)
plt.plot(e)
plt.title('Señal vocal "e"')
#Podemos apreciar la cuasiperiocidad de la señal
#El período de pitch lo calculamos directamente en la gráfica, sabiendo que disponemos de 256 muestras a 
#una fercuencia de muestreo de 8 kHz

#y el período de pitch 0.00621 s con una frecuencia de pitch de 161.03 Hz

#Ahora calculamos la autocorrelación
rx = np.correlate(e,e,'full') / N
k_=np.linspace(-N+1,N-1,2*N-1)
plt.figure(2)
plt.plot(k_,rx) #señal cuasi-periódica = autocorrelaación cuasi-periódica
plt.title('Autocorrelación')
#En la autocorrelación también podemos obtener el período de pitch

#Obtenemos el periodograma como la transformada de fourier de la autocorrelación
Px = np.abs(np.fft.fft(rx,Lfft))
w = np.linspace(-np.pi,np.pi,Lfft) #representamos de 0 a pi
plt.figure(3)
plt.semilogy(w,Px) #representamos en escala logaritmica
plt.grid()
plt.title('Periodograma')
#Podemos calcular también la frecuencia de pitch en el espectro como la diferencia
#entre dos picos (cualesquiera) sabiendo que el punto 1024 corresponde a la frecuencia de muestreo 8 kHz
# y vemos que coincide con la anteriormente calculada
#Si nos fijamos en la envolvente, podemos identificar también las frecuencias formantes
#que corresponden a los valores máximos o picos de la envolvente, son las frecuencias de
#resonancia del tracto vocal

#La frecuecia de pitch es 161.7 Hz y el período de pitch por tanto es de 0.00621
#Midiendo de forma aproximada las formantes, de 0 a pi observamos 3 resonancias, una en 445.63 Hz,
#otra en 1.757 kHz y otra en 3.514 kHz aproximadamente


#%% ESPECTRO LPC/AR

p = 12 #número de polos para el filtro
a = np.empty(12)

a_,sig2,k = aryule(e,p)
a = np.concatenate(([1],a_), axis = 0)
#Calculamos el espectro AR con la fórmula teórica
Px_AR = sig2 / (np.abs(np.fft.fft(a,Lfft)) ** 2)
plt.figure(4)
plt.semilogy(w,Px_AR)
plt.grid()
plt.title('PSD modelo AR')
plt.figure(5)
plt.semilogy(w,Px)
plt.semilogy(w,Px_AR)
plt.grid()
plt.title('Comparación PSDs')
#Al superponer las dos gráficas, nos damos cuenta que al representar el espectro
#AR con los coeficientes, nos proporciona la envolvente del periodograma
#Aquí se pueden apreciar mejor las frecuencias resonantes, las cuales
#coinciden con las frecuencias formantes

#Las frecuencias formantes, medidas en la envolvente son 445.63 Hz, 1.83 kHz y 3.501 kHz
#Comprobamos que coinciden con las calculadas aproximadamente únicamente con el periodograma

#%% ANÁLISIS HOMOMÓRFICO

#Conociendo la analogía entre potencias de salida y de entrada de un proceso, podemos
#introducir logaritmos, lo que se traduce en descomponer factores en sumandos, esto nos
#proporciona un valor llamado cepstrum, el cual diferencia la contribución de la excitación
#y la contribución del filtro log(Px) = log(|H|^2) + log(Pu)

cepstrum = np.fft.ifft(np.log(Px)) #Hacemos la transformada inversa para verlo en el dominio temporal

cepstrum_rep = cepstrum[1:100] #Representamos los 100 primeros coeficientes excluyendo el cx(0)
#Hacemos el módulo del cepstrum
c = np.empty(len(cepstrum_rep))
for i in range (len(cepstrum_rep)):
    c[i] = abs(cepstrum_rep[i])
    
    
plt.figure(7)

plt.plot(c)
plt.title('Cepstrum')
#La contribución del filtro la encontramos a valores bajos de cuefrencia y la contribución
#de la excitación la distinguimos en valores altos de cuefrencia. Las dos contribuciones, se 
#separan en el coeficiente 20 aproximadamente, el período de pitch vendrá dado por el valor máximo de 
#cuefrencia a partir de este coeficiente 20 (excitación). La entonación o el sexo del hablante estarán
#implícitos en la señal de excitación

#Medimos la frecuencia de pitch como Fs / n, nos da un valor de 167.01 Hz y el período de pitch
#es de 0.00598 s aproximadamente; vemos que coinciden con los valores del pitch anteriormente medidos

#Ahora vamos a obtener el cepstrum a partir del espectro LPC, el cual sólo tendrá información
#sobre el tracto vocal, y no sobre la excitación
cepstrum_LPC = np.fft.ifft(np.log(Px_AR))
cepstrum_LPC_rep = cepstrum_LPC[1:100]
c_LPC=np.empty(len(cepstrum_LPC_rep))

for i in range (len(cepstrum_LPC_rep)):
    c_LPC[i]=abs(cepstrum_LPC_rep[i])
    

plt.figure(8)
#plt.plot(cepstrum_LPC_rep)
plt.plot(c_LPC)
plt.title('Cepstrum LPC')

'''

plt.figure(9)
plt.plot(cepstrum_rep)
plt.plot(cepstrum_LPC_rep)
plt.title('Comparación cepstrum')
'''

plt.figure(9)
plt.plot(c)
plt.plot(c_LPC)
plt.title('Comparación cepstrum modulo')

#Como apreciamos en la figura anterior, solo es representada la señal que modela el 
#tracto vocal (filtro), ya que con el proceso AR hemos modelado la parte no sonora
#Vemos como, aproximadamente, en el coeficiente 20 (cuando comienza la señal de excitación)
#el cepstrum LPC pasa a ser nulo

#%% RECONOCIMIENTO DE VOCALES

#Calculamos el cepstrum de cada vocal, mediante la función 'vocal_x', con el mismo procedimiento descrito
#anteriormente, incluido el cepstrum para la vocal x.
cepstrum_a = vocal_x(1)
cepstrum_e = vocal_x(2)
cepstrum_i = vocal_x(3)
cepstrum_o = vocal_x(4)
cepstrum_u = vocal_x(5)
cepstrum_x = vocal_x(6)

#Hacemos una comparativa de todos los cepstrum
plt.figure(10)
plt.plot(cepstrum_a, label = 'vocal a')
plt.plot(cepstrum_e, label = 'vocal e')
plt.plot(cepstrum_i, label = 'vocal i')
plt.plot(cepstrum_o, label = 'vocal o')
plt.plot(cepstrum_u, label = 'vocal u')
plt.plot(cepstrum_x, label = 'vocal x')
plt.legend()
plt.title('Cepstrum vocales')
#Mediante la función 'indice_vocal' calculamos la distancia mínima
indice_vocal = distancia_minima(cepstrum_a, cepstrum_e, cepstrum_i, cepstrum_o, cepstrum_u, cepstrum_x)
#Ahora imprimimos el resultado en pantalla
identifica_vocal(indice_vocal)
print(indice_vocal)
    

