
clc;
clear all;
close all;

%Determinamos variables
Fs = 16e3;    %Frecuencia de muestreo
d = 0.04;     %Distancia entre sensores
Vprop = 340;  %Velocidad del sonido
N = 7;        %Contamos con 7 elementos
phi = pi/4;   %Ángulo de llegada del target
iter = 1;        %Iterador del bucle para las ventanas
n = (0:1:6);         
tn = ((d*cos(phi).*n)/Vprop);%Creamos el vector de retardos
theta = linspace(0, pi, 200); % Barrido en theta
theta_polar = linspace(0, 2*pi, 400); % Barrido representación polar

w = pesos(tn);

D = calcula_Dteo(w, 100, d, Vprop,theta);

figure
polarplot(theta_polar, D);