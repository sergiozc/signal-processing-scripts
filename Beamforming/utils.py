# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 21:36:09 2023

@author: Usuario
"""

#%% PATRÓN DE DIRECTIVIDAD EXPERIMENTAL


A = np.zeros(len(barrido))
barrido_exp = np.linspace(np.pi / 2, 3 * np.pi / 2, 36)
barrido_exp_sim = np.linspace(np.pi / 2, 5 * np.pi / 2, 72)
barrido_real = np.linspace(0, np.pi, 36)

pos = 0
for i in range(len(barrido)):
    valores = np.abs(y[pos:(pos + segmento)])
    valores = np.sort(valores)
    valores = valores[::-1]
    A[i] = np.mean(valores[:30000])
    pos = pos + segmento

Amax = max(A)

#Se da la vuelta a las amplitudes por hacer el barrido de 270 hasta 90
Dexp = A[::-1] / Amax

#Parte simétrica de la directividad (se ha hecho la prueba desde 90 grados)
Dsim_exp = Dexp[::-1]
Dtot_exp = np.concatenate((Dexp, Dsim_exp), axis=0)


plt.figure(18)
plt.polar(barrido_exp, abs(Dexp))
plt.title('Diagrama directividad experimental')

plt.figure(19)
plt.polar(barrido_exp_sim, abs(Dtot_exp))
plt.title('Diagrama directividad experimental')
    



#%% MSE directividad

suma = 0
for phi_ in range (len(barrido_polar)):
    suma = suma + (Dtot[phi_] - Dtot_exp[phi_])

MSE = suma **2 / len(barrido_polar)

print('El MSE de la directividad es: ' + str(MSE))

#%%Factor de directividad
denom = 0
num = (abs(Dtot_exp[0]))**2

for phi2 in range (len(barrido_polar)):
    denom = (abs(Dtot_exp[phi2]))**2 + denom
    
denom = denom / len(barrido_polar)
Fd = (num / denom)

#%% Representación

plt.figure(20)
plt.polar(barrido_polar, abs(Dtot))
plt.polar(barrido_exp_sim, abs(Dtot_exp))
plt.title('Comparación')