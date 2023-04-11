
clc;
clear all;
close all;


% BEAMFORMING DELAY AND SUM
%Determinamos variables
Fs = 16e3;    %Frecuencia de muestreo
d = 0.04;     %Distancia entre sensores
Vprop = 340;  %Velocidad del sonido
N = 7;        %Contamos con 7 elementos
phi = pi/4;   %Ángulo de llegada del target
n = (0:1:6);  %Vector que computa el núemero de canal       
tn = (d*cos(phi).*n)/Vprop;%Creamos el vector de retardos
theta = linspace(0, pi, 200); % Barrido en theta
theta_polar = linspace(0, 2*pi, 400); % Barrido representación polar
freq = linspace(1, 8000, 129);



w = pesos(tn, freq);

% Visualización de los retardos impuestos
% (La señal llega antes al último sensor).
figure(1)
stem(tn)
title('Retardos')
ylabel('tiempo (s)')
xlabel('N sensor')



%% Directividad para 1 frecuencia
index_freq1 = 17;
D_1kHz = calcula_Dteo(w, freq, index_freq1, d, Vprop, theta);
index_freq2 = 129;
D_8kHz = calcula_Dteo(w, freq, index_freq2, d, Vprop, theta);


figure(2);
polarplot(theta_polar, D_1kHz);
hold on
polarplot(theta_polar, D_8kHz);
title('Directividad teórica');
legend('f = 1kHz', 'f = 8kHz');

%% Directividad en función de la frecuencia
theta_surf = linspace(0, 2*pi, 129); % Barrido en theta
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


