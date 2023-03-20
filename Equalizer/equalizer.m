%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PRÁCTICA 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Ecualizadores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Javier Lobato Martín
% Sergio Zapata Caparrós


%% Ecualización Fija

clc;
close all;
clear all;

% Hacemos un bucle que nos vaya comparando los distintos valores de u para
% encontrar el óptimo

N = 300;        % Tamaño de la secuencia de entrenamiento (filas de la matriz)
%u = 2;          % Retardo de ecualización (determinar el óptimo)
offset = 128;   % Rango 0-255
M = 8;          % Número de coeficientes del filtro (columnas de la matriz)
SNR_array = zeros(1, M+1); % Vector que determinará la mejor SNR según u
SER_array = zeros(1, M+1); % Vector que determinará la mejor SER según u

xi = load('Tx_fijo_variables.mat');

s = xi.trainseq; % Secuencia de entrenamiento
z = xi.z;       % Señal de salida del canal total (canal físico + filtro receptor)
ini_z = 101;     % Donde comienza la señal z (después de los 100 símbolos aleatorios)
origi = xi.x;   % Señal de audio original
%soundsc(origi);

scatterplot(z);
title('Señal de salida canal total');

for u = 0:8

% Primero creamos la matriz Z en función del retardo u
Z_matrix = zeros(N, M);

for i = 1:N
    for j = 1:M+1
        Z_matrix(i, j) = z(ini_z - j + i + u);
    end  
end

% Con la matriz, hayamos los coeficientes del filtro (conocemos la expresión óptima del filtro)
f = inv(Z_matrix' * Z_matrix) * (Z_matrix' * s);

% La señal deseada
d = static_filter(f, z);
%d = filter(f, 1, z); % Filtrado fijo con la función filter()
scatterplot(d);
title("Señal de salida ecualizador retardo u = " + u);

% Hay que tener en cuenta el retardo que introduce el parámetro u
d = circshift(d,-u);

% Demodulación y decodificación
demodulada = qamdemod(d, 16); % Demodulador 16-QAM
xsint = qam16_to_audio8bit(demodulada); % Separación de MSB y LSB
%soundsc(xsint); % Señal resintetizada

% Se eliminan las 100 muestras añadidas de relleno (200 después de la
% decodificación por LSB y MSB) y se hace un vector columna
xsint = transpose(xsint(201:end));

% Cálculo de la SNR
pot_signal = sum(origi.^2);
pot_ruido = sum((xsint - origi).^2);
SNR = calculo_SNR(pot_signal, pot_ruido);
fprintf('SNR(TSE-LS)  = %f dB para u = %i \n', SNR, u);


% Cálculo del Symbol Error Rate (SER)
% demodulada contiene la secuencia de símbolos en formato 0-16
% origi contiene la señal original de audio, debemos pasarla al formato de
% símbolos para la comparación

origi_symb = audio8bit_to_qam16(origi); %codificamos QAM el audio original
demodulada = demodulada(401:end); %igualamos las longitudes

% comparación es un vector binario que indica si los elementos son iguales
% o no entre las dos secuencias
comparacion = demodulada-origi_symb'; 
comparacion = comparacion./comparacion;
comparacion(isnan(comparacion)) = 0;

SER = sum(comparacion)/length(comparacion)*100; %calculamos el SER
fprintf('SER(TSE-LS) (%%) = %f %% para u = %i \n', SER, u);

% Muestra de resultados para u
SNR_array(u+1) = SNR;
SER_array(u+1) = SER;

end

% soundsc(xsint);

u_array = (0:u) ;
figure
stem(u_array, SNR_array)
title('SNR frente al retardo u')
xlabel('parámetro u')
ylabel('SNR (dB)')

figure
stem(u_array, SER_array)
title('SER frente al retardo u')
xlabel('parámetro u')
ylabel('SER')
ylim([0 12]);

%Como vemos, el parámetro u ideal es 2. Esto tiene sentido ya que coincide
%exactamente con L/2, siendo L la longitud de un símbolo (4 muestras)

%% Ecualización adaptable
clc;
close all;
clear all;

%mu_vector = [0.0025 0.005, 0.01, 0.02, 0.03]; %Creamos un vector para probar distintos valores del paso en el descenso de gradiente
mu_vector = linspace(0.0025,0.03,5); %Creamos un vector para probar distintos valores del paso en el descenso de gradiente
SNR_vector = zeros(1, length(mu_vector));
SER_vector = zeros(1, length(mu_vector));


for k = 1:length(mu_vector)

    N = 300;        % Tamaño de la secuencia de entrenamiento (filas de la matriz)
    u = 2;          % Retardo de ecualización
    offset = 128;   % Rango 0-255
    M = 8;              % Número de coeficientes del filtro
    %mu = 0.01;          % Paso correspondiente al avance en descenso de gradiente
    mu = mu_vector(k);

    xi = load('Tx_var_variables.mat');
    z = xi.z;           % Señal de salida del canal total (canal físico + filtro receptor)
    s = xi.trainseq;
    origi = xi.x;
    ini_z = 101;     % Donde comienza la señal z (después de los 100 símbolos aleatorios)
    
    
    
    % Primero creamos la matriz Z en función del retardo u
    Z_matrix = zeros(N, M);
    
    for i = 1:N
        for j = 1:M+1
            Z_matrix(i, j) = z(ini_z - j + i + u);
        end  
    end
    
    % Con la matriz, hayamos los coeficientes del filtro (conocemos la expresión óptima del filtro)
    f = inv(Z_matrix' * Z_matrix) * (Z_matrix' * s);
    
    
    d = adaptable_filter(f, z, mu); % Filtrado adaptable (ver función)
    
    scatterplot(d);
    title("Señal de salida ecualizador retardo u = " + u + " y mu = " + mu);
    
    
    % Hay que tener en cuenta el retardo que introduce el parámetro u
    d = circshift(d,-u);
    
    % Demodulación y decodificación
    demodulada = qamdemod(d, 16); % Demodulador 16-QAM
    xsint = qam16_to_audio8bit(demodulada); % Separación de MSB y LSB
    soundsc(xsint); % Señal resintetizada
    
    % Se eliminan las 100 muestras añadidas de relleno (200 después de la
    % decodificación por LSB y MSB) y se hace un vector columna
    xsint = transpose(xsint(201:end));
    
    % Cálculo de la SNR
    pot_signal = sum(origi.^2);
    pot_ruido = sum((xsint - origi).^2);
    SNR = calculo_SNR(pot_signal, pot_ruido);
    fprintf('SNR(NMLS)  = %f dB para u = %i y mu = %i \n', SNR, u, mu);
    
    
    % Cálculo del Symbol Error Rate (SER)
    % demodulada contiene la secuencia de símbolos en formato 0-16
    % origi contiene la señal original de audio, debemos pasarla al formato de
    % símbolos para la comparación
    
    origi_symb = audio8bit_to_qam16(origi); %codificamos QAM el audio original
    demodulada = demodulada(401:end); %igualamos las longitudes
    
    % comparación es un vector binario que indica si los elementos son iguales
    % o no entre las dos secuencias
    comparacion = demodulada-origi_symb'; 
    comparacion = comparacion./comparacion;
    comparacion(isnan(comparacion)) = 0;
    
    SER = sum(comparacion)/length(comparacion)*100; %calculamos el SER
    fprintf('SER(NMLS) (%%) = %f %% para u = %i y mu = %i \n', SER, u, mu);
    

    
    SNR_vector(k) = SNR;
    SER_vector(k) = SER;
end

%Representamos la evolución de nuestros indicadores
figure
stem(mu_vector,SNR_vector)
title('Evolución de la SNR frente al parámetro \mu')
xlabel('\mu')
ylabel('SNR (dB)')

figure
stem(mu_vector,SER_vector)
title('Evolución del SER frente al parámetro \mu')
xlabel('\mu')
ylabel('SER')
%Como podemos comprobar, el parámetro mu que optimiza ambas medidas (SNR y SER)
%se encuentra en torno a mu = 0.01
%Es un paso mucho más cercano a 0 que al límite superior del parámetro de
%convergencia.

%El retardo u de ecualización va a tener el mismo valor que en la primera
%sección, ya que el valor óptimo es L/2, siendo L el tamaño del símbolo,
%que es 4 muestras en este caso.


