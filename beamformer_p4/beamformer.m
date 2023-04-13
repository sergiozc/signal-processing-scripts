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
Ntramas = L_signal/128;       %Determinamos el número de tramas
iter = 1;                     %Iterador del bucle para las ventanas
win = hann(Ltrama+1);         %Establecemos la ventana de Hanning
n = (0:1:6);         
tn = ((d*cos(phi).*n)/Vprop);  %Creamos el vector de retardos
freq = linspace(1, 8000, 129); %Vector de frecuencias


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
        Xn = Xn(1:129);                 %Tomamos las componentes de frecuencia de 0 a Fs/2 (Fs/2 = 8 kHz)s     
        Xn = Xn .* w(:,c);              %Multiplicamos por los pesos correspondientes
        

        
        
        
        % ESTO CREO QUE VA FUERA DEL FOR (DESPUÉS DEL SUM)
        %Realizamos la simetrización para practicar la transformada inversa
        XOUT = cat(1, Xn, Xn(end:-1:2)); %NO SE SI ESTÁ BIEN LA DIMENSIÓN

        %Hacemos la transformada de Fourier inversa
        xout = ifft(sqrt(win).*XOUT);
        xout = real(xout);
        
        %Sacamos el resultado final (concatenación de las tramas)
        xc_out(iter:iter + Ltrama, c) = xc_out(iter:iter + Ltrama, c) + xout;

    end
    
    
    iter = iter + 127; %Actualizamos el iterador

end


%% comprobaciones

figure
plot(xc_out(:,3))

xc_out_all = zeros(L_signal,1);

%Sumando todos (SUM)

for i =1:7

    xc_out_all = xc_out_all + xc_out(:,i);

end

figure
plot(xc_out_all)
title('suma de todos ')


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

