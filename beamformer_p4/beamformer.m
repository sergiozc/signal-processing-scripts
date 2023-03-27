%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PRÁCTICA 5: BEAMFORMER DELAY AND SUM %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autores: Javier Lobato Martín y Sergio Zapata Caparrós
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

%Cargamos las señales
signals = load('signals_array.mat');
xc = signals.xc;
xc = cell2mat(xc); % 7 columnas

%Determinamos variables
Fs = 16e3;    %Frecuencia de muestreo
d = 0.04;     %Distancia entre sensores
Vprop = 340;  %Velocidad del sonido
Ltrama = 256; %Tramas de 256 muestras
N = 7;        %Contamos con 7 elementos
phi = pi/4;   %Ángulo de llegada del target
L_signal = length(xc(:,1)); %Longitud total de la señal
Ntramas = L_signal/128;     %Determinamos el número de tramas
iter = 1;        %Iterador del bucle para las ventanas
win = hann(Ltrama+1);         %Establecemos la ventana de Hanning
n = [0:1:6];         
tn = ((d*cos(phi).*n)/Vprop);%Creamos el vector de retardos

pesoss = pesos(tn);


%Para realizar el trabajo, usaremos un doble bucle
for ntram = 1:Ntramas

    

    for c = 1:7
        
        
        xn = xc(iter:iter + Ltrama ,c); %Tomamos la porción de señal del canal correspondiente
        Xn = fft(sqrt(win).*xn); %Realizamos la transformada de Fourier de la ventana
        Xn = Xn(1:129);          %Tomamos las componentes de frecuencia de 0 a Fs/2 (Fs/2 = 8 kHz)
        




        %Practicamos la FFT
        
        




    end
    
    iter = iter + 127; %Actualizamos el iterador

end

