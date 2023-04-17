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
Fs = 16e3;          %Frecuencia de muestreo
d = 0.04;           %Distancia entre sensores
Vprop = 340;        %Velocidad del sonido
Ltrama = 256;       %Tramas de 256 muestras
N = 7;              %Contamos con 7 elementos
phi = pi/4;         %Ángulo de llegada del target
L_signal = length(xc(:,1));   %Longitud total de la señal
Ntramas = floor(L_signal/128);       %Determinamos el número de tramas
iter = 1;                     %Iterador del bucle para las ventanas
win = hanning(Ltrama+1,'periodic');
n = (0:1:6);         
tn = ((d*cos(phi).*n)/Vprop);  %Creamos el vector de retardos
freq = linspace(1, 8000, 129); %Vector de frecuencias (Fs >= Fmax)


%Creamos la matriz vacía donde vamos a guardar el resultado final
xc_out = zeros(L_signal,N);

%obtenemos los pesos por medio de nuestra función auxiliar
w = pesos(tn, freq);
XOUT = zeros(129, 1);

%Para realizar el trabajo, usaremos un doble bucle
for ntram = 1:Ntramas 

    for c = 1:7        
        
        xn = xc(iter:iter + Ltrama ,c); %Tomamos la porción de señal del canal correspondiente
        Xn = fft(sqrt(win).*xn);        %Realizamos la transformada de Fourier de la ventana
        Xn = Xn(1:Ltrama/2+1);          %Tomamos las componentes de frecuencia de 0 a Fs/2 (Fs/2 = 8 kHz)s     
        Xn = Xn .* conj(w(:,c));        %Multiplicamos por los pesos correspondientes
        

        
        %Realizamos la simetrización para practicar la transformada inversa
        XOUT = cat(1, Xn, conj(Xn(end:-1:2)));
        xout = real(ifft(sqrt(win).*XOUT));
        
        %Concatenación de tramas mediante ''overlap add''
        xc_out(iter:iter + Ltrama, c) = xc_out(iter:iter + Ltrama, c) + xout;

    end
    
    
    iter = iter + 127; %Actualizamos el iterador

end

xc_out_sum = sum(xc_out, 2);
soundsc(real(xc_out_sum),Fs);


figure(5)
plot(real(xc_out_sum));


%% Cálculo SNR

% SNR DESPUÉS DEL BEAMFORMING DAS
ruido_orig = var((xc(1:3000, 1))); %Interferencia aislada en las 3000 primeras muestras
pot_orig = var((xc(3001:end, 1)));
SNR_orig = calculo_SNR(pot_orig, ruido_orig);
fprintf('SNR(antes)  = %f dB\n', SNR_orig);

% SNR DESPUÉS DEL BEAMFORMING DAS
ruido_DAS = var(real(xc_out_sum(1:3000)));
pot_DAS = var(real(xc_out_sum(3001:end)));
SNR_DAS = calculo_SNR(pot_DAS, ruido_DAS);
fprintf('SNR(después)  = %f dB\n', SNR_DAS);


%% Comprobación beamformer

% VARIABLES PARA REPRESENTAR
theta = linspace(0, pi, 200); % Barrido en theta
theta_polar = linspace(0, 2*pi, 400); % Barrido representación polar
theta_surf = linspace(0, 2*pi, 129); % Barrido en theta



% VISUALIZACIÓN DE LOS RETARDOS
figure(1)
stem(tn)
title('Retardos')
ylabel('tiempo (s)')
xlabel('N sensor')
% (La señal llega antes al último sensor).


% DIRECTIVIDAD PARA 1 FRECUENCIA
index_freq1 = 17; % 1kHz
D_1kHz = calcula_Dteo(w, freq, index_freq1, d, Vprop, theta);
index_freq2 = 129; % 8kHz
D_8kHz = calcula_Dteo(w, freq, index_freq2, d, Vprop, theta);

figure(2);
polarplot(theta_polar, D_1kHz);
hold on
polarplot(theta_polar, D_8kHz);
title('Directividad teórica');
legend('f = 1kHz', 'f = 8kHz');

% DIRECTIVIDAD EN FUNCIÓN DE LA FRECUENCIA
Df = calcula_Df(w, freq, d, Vprop, theta_surf);
figure(3);
surf(freq, rad2deg(theta_surf), Df);
xlabel('f(Hz)');
ylabel('phi');
zlabel('D(f, phi)');

figure(4);
pcolor(freq, rad2deg(theta_surf), abs(Df.'));
shading interp;
colorbar;
xlabel('Frecuencia (Hz)');
ylabel('φ(grados)');
title('Directividad en función de la frecuencia y ángulo');
% Se aprecia como a frecuencias más altas, el ancho del lóbulo principal
% disminuye (se hace más directivo). delta = 1/(f*L)
% Se ve también que a partir de un determinado valor de frecuencia,
% aparecen más rizados (o lóbulos secundarios).

