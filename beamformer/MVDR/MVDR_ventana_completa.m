clc;
clear all;
close all;


%% Lectura de las señales

% DATOS
% Array lineal de espaciado variable 4cm (4A)
% Posicion Speaker: 1 m del centro del array
% Muestras: 16 kHz, 16 bits por muestra
% 15 canales
Fs     = 16000; % Frec. muestreo
Narray = 15; % Nº de canales del array
dist=[0 16 24 32 36 40 44 48 52 56 60 64 72 80 96]*0.01; % Espaciado (m)
c=340; % Velocidad propagacion
fm = 16e3;

% Seleccionar señales
dir = 'signals/';
fname = 'an101-mtms-arr4A.adc'; dspk=1; % computer lab, speaker 1m

% Lectura de las señales
fnamebase=fname(1:16);
fname = strcat(dir,fname);
[fid,msg] = fopen(fname,'r','b');
if fid < 0
  disp(msg);
  exit;
else
  data = fread(fid,'int16');
  fclose(fid);
end

% Separa señales
nsamp=[]; x={};
for i = 1:Narray
    x{i} = data(i:Narray:end);
    x{i} = x{i} - mean(x{i});
    nsamp(i)=length(x{i});
end

% Seleccionamos subarray
index=[5, 6, 7, 8, 9, 10, 11]; %array de 4 cm
%index=[3, 4, 6, 8, 10, 12, 13]; %array de 8 cm
%index=[1, 2, 4, 8, 12, 14, 15]; %array de 16 cm
%index=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]; %array completo
Nc=length(index); % No. de canales a utilizar
dist=dist-dist(index(1)); %Primer elemento subarray como referencia
dist=dist(index);  %fijar espaciado subarray

% Agrupamos señales subarray en matriz
Nsamp=min(nsamp); % No. de muestras
y=[];
for n=1:Nc
    y=[y x{index(n)}(1:Nsamp)];
end
maxmax=max(max(abs(y)));
y=y/maxmax; %normalización del rango de la señal

% Reproduce y guarda central como señal de referencia ruidosa
ncent=floor(Nc/2)+1;
xcent=y(:,ncent);
%sound(xcent,fm)
fcent=strcat(fnamebase,'.wav');
audiowrite(fcent,xcent/max(abs(xcent)),Fs)


figure(1)
plot(y(:,ncent))
title('Representación temporal del audio')

%% Definición de variables

Fs = 16e3;          % Frecuencia de muestreo
d=dist(2)-dist(1);  % Separación entre elementos del array
Vprop = 340;        % Velocidad del sonido
Ltrama = 256;       % Tramas de 256 muestras
Lfft = 512;         % Longitud de la FFT
N = Nc;             % Contamos con 7 elementos
phi = pi/2;         % Ángulo de llegada del target
L_signal = length(y(:,ncent));   %Longitud total de la señal
win = hanning(Ltrama+1,'periodic'); % Ventana de hanning  
freq = linspace(0,256,257)*(Fs/Lfft); % de 0 a 8000 Hz
n=0:1:N-1; % Índice de los elementos del array
c = 340; % Velocidad de propagación

%% Tipo de onda
% Se puede elegir entre onda plana o esférica
[d_n, tn] = onda_tipo(c, d, N, n, 'spherical');

%% MATRIZ DE CORRELACIÓN ESPACIAL DEL RUIDO
muestras_ruido = 8000;
corr_noise = noise_matrix_vent_complet(N, freq, win, Ltrama, Lfft, muestras_ruido, y);

%% Cómputo del Beamformer

% CÁLCULO DE PESOS
% Pesos MVDR
w = pesos_MVDR(d_n, tn, freq, corr_noise);fprintf('Beamformer: MVDR \n');
% Pesos DAS
%w = pesos_DAS(d_n,tn, freq);fprintf('Beamformer: DAS \n');

% SEÑAL DIVISIBLE ENTRE Ltrama 
[m,~]=size(y);
resto=mod(m,Ltrama);
y=y(1:m-resto,:);
% Se obtiene el número de muestras que tendrá la señal sobre la que se
% aplicará el beamforming
[m,~]=size(y); 

Ntramas=2*(m/Ltrama)-1;

%% PROCESADO POR TRAMAS (análisis-síntesis)

xc_out = zeros(L_signal,N); % Matriz del resultado final
XOUT = zeros(Lfft/2+1, 1); % Señal depúes de aplicar los pesos
iter = 1;
for ntram = 1:Ntramas  % Se computa cada trama

    for c = 1:N        % Se computa cada sensor
        
        xn = y(iter:iter + Ltrama ,c); %Tomamos la porción de señal del canal correspondiente
        Xn = fft(win.*xn, Lfft);        %Realizamos la transformada de Fourier de la ventana (512 muestras)
        Xn = Xn(1:Lfft/2+1);          %Tomamos las componentes de frecuencia de 0 a Fs/2 (Fs/2 = 8 kHz)s     
        Xn = Xn .* conj(w(:,c));        %Multiplicamos por los pesos correspondientes
        

        %Realizamos la simetrización para practicar la transformada inversa
        simet = conj(Xn);
        XOUT = cat(1, Xn, simet(end:-1:2));
        xout = real(ifft(XOUT));
        
        %Concatenación de tramas mediante ''overlap add''
        xc_out(iter:iter + Lfft, c) = xc_out(iter:iter + Lfft, c) + xout;

    end
    
    iter = iter + 127;
end

%Unimos todos los canales y realizamos una escucha
xc_out_sum = sum(xc_out, 2);
% Eliminamos la cola residual de la ultima trama
xc_out_sum=xc_out_sum(1:end-Lfft/2);

% Normalizamos la señal y la escuchamos
xout_norm = xc_out_sum/max(abs(xc_out_sum));
soundsc(real(xout_norm),Fs);

% Guardamos señal resultante normalizada
fout=strcat('Resultado','.wav');
audiowrite(fout,xout_norm,Fs)


figure(2)
plot(y(:,ncent));
hold on
plot(real(xout_norm));
hold off
legend('Señal sensor central','Señal a la salida del beamformer')
title('Representación temporal tras MVDR')
%Se puede comprobar como el ruido se ha minimizado

%% Cálculo SNR
%Para realizar el cálculo de la SNR, calculamos la potencia de la señal
%y del ruido (primeras 8000 muestras) y obtenemos el ratio.

% SNR ANTES DEL BEAMFORMING
ruido_orig = var((y(1:8000, 1))); %Interferencia aislada en las 3000 primeras muestras
pot_orig = var((y(8000:end, 1)));
SNR_orig = calculo_SNR(pot_orig, ruido_orig);
fprintf('SNR(before)  = %f dB\n', SNR_orig);

% SNR DESPUÉS DEL BEAMFORMING
ruido_BF = var(real(xout_norm(1:8000)));
pot_BF = var(real(xout_norm(8000:end)));
SNR_BF = calculo_SNR(pot_BF, ruido_BF);
fprintf('SNR(after)  = %f dB\n', SNR_BF);


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
polarplot(theta_polar, D_matrix(3,:),'Color',[color_map(3,:)]);
polarplot(theta_polar, D_matrix(4,:),'Color',[color_map(4,:)]);
polarplot(theta_polar, D_matrix(5,:),'Color',[color_map(5,:)]);
polarplot(theta_polar, D_matrix(6,:),'Color',[color_map(6,:)]);
polarplot(theta_polar, D_matrix(7,:),'Color',[0.6350 0.0780 0.1840]);
polarplot(theta_polar, D_matrix(8,:),'Color',[color_map(7,:)]);
polarplot(theta_polar, D_matrix(9,:),'Color',[color_map(8,:)]);
polarplot(theta_polar, D_matrix(10,:),'Color',[color_map(9,:)]);
polarplot(theta_polar, D_matrix(11,:),'Color',[color_map(10,:)]);
hold off
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

