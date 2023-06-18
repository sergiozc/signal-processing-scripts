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

figure(1)
plot(xc(:,3))
title('Representación temporal del audio')


%Determinamos variables
Fs = 16e3;          %Frecuencia de muestreo
d = 0.04;           %Distancia entre sensores
Vprop = 340;        %Velocidad del sonido
Ltrama = 256;       %Tramas de 256 muestras
Lfft = 512;
N = 7;              %Contamos con 7 elementos
phi = pi/4;         %Ángulo de llegada del target
L_signal = length(xc(:,1));   %Longitud total de la señal
Ntramas = floor(L_signal/128);       %Determinamos el número de tramas
iter = 1;                      %Iterador del bucle para las ventanas
win = hanning(Ltrama+1,'periodic'); %Establecemos la ventana de hanning
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

    for c = 1:N        
        
        xn = xc(iter:iter + Ltrama ,c); %Tomamos la porción de señal del canal correspondiente
        Xn = fft(win.*xn);        %Realizamos la transformada de Fourier de la ventana
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

%Unimos todos los canales y realizamos una escucha
xc_out_sum = sum(xc_out, 2);
soundsc(real(xc_out_sum),Fs);

xout_norm = xc_out_sum/max(abs(xc_out_sum));
x_antes_norm = xc(:,3)/max(abs(xc(:,3)));

% Guardamos señal resultante normalizada
fout=strcat('Resultado','.wav');
audiowrite(fout,xout_norm,Fs)


figure(2)
plot(x_antes_norm);
hold on
plot(real(xout_norm));
hold off
legend('Señal sensor central','Señal a la salida del beamformer')
title('Representación temporal tras D&S')
%Se puede comprobar como el ruido se ha minimizado

%% Cálculo SNR
%Para realizar el cálculo de la SNR, calculamos la potencia de la señal
%y del ruido (primeras 3000 muestras) y obtenemos el ratio.

% SNR DESPUÉS DEL BEAMFORMING D&S
ruido_orig = var((xc(1:3000, 1))); %Interferencia aislada en las 3000 primeras muestras
pot_orig = var((xc(3001:end, 1)));
SNR_orig = calculo_SNR(pot_orig, ruido_orig);
fprintf('SNR(antes)  = %f dB\n', SNR_orig);

% SNR DESPUÉS DEL BEAMFORMING DAS
ruido_DAS = var(real(xc_out_sum(1:3000)));
pot_DAS = var(real(xc_out_sum(3001:end)));
SNR_DAS = calculo_SNR(pot_DAS, ruido_DAS);
fprintf('SNR(D&S)  = %f dB\n', SNR_DAS);


%% Comprobación beamformer

% VARIABLES PARA REPRESENTAR
theta = linspace(0, pi, 200); % Barrido en theta
theta_polar = linspace(0, 2*pi, 400); % Barrido representación polar
theta_surf = linspace(0, 2*pi, 129); % Barrido en theta


% VISUALIZACIÓN DE LOS RETARDOS
figure(3)
stem(tn)
title('Retardos de cada sensor')
ylabel('tiempo (s)')
xlabel('N sensor')
% La señal llega primero al último sensor (7). Es por ello que tiene
% asociado el mayor retardo de todos.


% CÁLCULO DE DIRECTIVIDADES
% Almacenamos en un vector los índices correspondientes a las frecuencias
% de interés 

% frecuencias = [100, 400, 700, 1000, 2000,
% 3000, 4000, 5000, 6000, 7000, 8000]
index_freq = [3, 7, 12, 17, 33, 49, 69, 81, 97, 113, 129];
D_matrix = zeros(length(index_freq),400);

for index=1:length(index_freq)

    D_matrix(index,:) = calcula_Dteo(w, freq, index_freq(index), d, Vprop, theta);

end

index_freq1 = 17; % 1kHz
D_1kHz = calcula_Dteo(w, freq, index_freq1, d, Vprop, theta);
index_freq2 = 129; % 8kHz
D_8kHz = calcula_Dteo(w, freq, index_freq2, d, Vprop, theta);

% Crear un gradiente de 10 colores que va de verde a azul
color_map = [linspace(0,0,5)', linspace(1,0,5)', linspace(1,1,5)';
             linspace(0,0,5)', linspace(0,1,5)', linspace(1,0,5)'];

% Establecer el gradiente de colores personalizado
colormap(color_map);



figure(4);
polarplot(theta_polar, D_matrix(1,:),'Color',[color_map(1,:)]);
hold on
polarplot(theta_polar, D_matrix(2,:),'Color',[color_map(2,:)]);
hold on
polarplot(theta_polar, D_matrix(3,:),'Color',[color_map(3,:)]);
hold on
polarplot(theta_polar, D_matrix(4,:),'Color',[color_map(4,:)]);
hold on
polarplot(theta_polar, D_matrix(5,:),'Color',[color_map(5,:)]);
hold on
polarplot(theta_polar, D_matrix(6,:),'Color',[color_map(6,:)]);
hold on
polarplot(theta_polar, D_matrix(7,:),'Color',[0.6350 0.0780 0.1840]);
hold on
polarplot(theta_polar, D_matrix(8,:),'Color',[color_map(7,:)]);
hold on
polarplot(theta_polar, D_matrix(9,:),'Color',[color_map(8,:)]);
hold on
polarplot(theta_polar, D_matrix(10,:),'Color',[color_map(9,:)]);
hold on
polarplot(theta_polar, D_matrix(11,:),'Color',[color_map(10,:)]);
title('Directividad teórica');
legend('f = 100Hz','f = 400Hz','f = 700Hz','f = 1kHz','f = 2kHz','f = 3kHz','f = 4.25kHz','f = 5kHz','f = 6kHz','f = 7kHz', 'f = 8kHz');
% Se observa el aumento de la directividad conforme la frecuencia va
% creciendo.
% Podemos ver como aparecen lóbulos indeseados en las frecuencias a partir
% de nuestra frecuencia óptima de lambda/2 = 4250 Hz. El valor óptimo es
% donde se puede observar una mayor directividad sin lóbulos secundarios
% indeseados. Esto se justifica en base al teorema de muestreo espacial:
% La separación entre sensores 


% DIRECTIVIDAD EN FUNCIÓN DE LA FRECUENCIA
Df = calcula_Df(w, freq, d, Vprop, theta_surf);
figure(5);
surf(rad2deg(theta_surf), freq, Df);
ylabel('f(Hz)');
xlabel('phi');
zlabel('D(f, phi)');
title('Directividad en función de la frecuencia y ángulo');


figure(6);
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

