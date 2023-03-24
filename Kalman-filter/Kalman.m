%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PRÁCTICA 4: Filtro de Kalman %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autores: Javier Lobato Martín y Sergio Zapata Caparrós
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definición de variables

clc;
clear all; 
close all;

dist = 35786;   % Distancia en km a la superficie terrestre
v = 11068;      % Velocidad orbital (km/h)

x = [42164 0 0 11068];         % Vector de estado [x y x(.) y(.)]
                               % siendo y(.) la velocidad ('y' posición)
Tpred = 0.1;    % Tiempo de predicción (0.1 h)
Tmed = 3;       % Tiempo de medidas (3 h)
var_pos = 50;   % Varianza de error posición (km2)
var_v = 1;      % Varianza de error velocidad (km2/h2)
global Q        % Matriz covarianza error
global GM;      % Cte gravitacional * Masa de la Tierra


variables = load('Variables_Orbita.mat');
sig2w = variables.sig2w;    % Varianza de la Perturbación de la aceleración
GM = variables.GM;

w = [0 0 sig2w sig2w];      % Vector de las perturbaciones
Q = diag(w);                % Matriz Covarianza(de las perturbaciones)

P = diag([var_pos var_pos var_v var_v]); % Matriz de errores inicial (diagonal)

R = 0.1; % Error de medida

%% Construcción del filtro

% DEFINICIÓN DE VARIABLES

%Primero de todo, conformamos el vector xp inicial, concatenación de x y los
%elementos de la matriz P reorganizados:
xp = [x reshape(P,[1,16])];
orbita = zeros(2,length(variables.tt));
med = 1; %Contador que nos va a indicar qué medida tenemos que tomar
k = 1;   %Contador que nos va a indicar en qué instante nos encontramos (en formato entero)
error_cometido = zeros(1,length(variables.tt)); %Almacenaremos el error en las 240 iteraciones
error_estimado = zeros(1,length(variables.tt)); %Almacenaremos el error estimado en las 240 iteraciones
pos = zeros(2,length(variables.tt));   %Almacenaremos la posición estimada en las 240 iteraciones
r_est = zeros(1,length(variables.tt)); %Almacenaremos el radio estimado de órbita
r_est_medidas = zeros(1,240); %Aquí almacenaremos las medidas tomadas cada 3 horas (para comparar luego)
r_orig = zeros(1, 240); %Aquí vamos a almacenar el radio verdadero


for kt = 0:Tpred:24-Tpred

% PREDICCIÓN
%Para la primera iteración, tenemos ya conformado el vector xp incial. Para
%el resto de iteraciones, tenemos el vector xp de salida de la iteración
%anterior

% Límites resolución ecuaciones diferenciales
tini = kt;
tfin = tini + Tpred;

%Resolvemos con ode45 las ecuaciones diferenciales
[~, xp] = ode45(@difeq, [tini, tfin], xp);


xp = xp(end, :); % Tomamos las 4 variables de estado


% INCORPORACIÓN DE MEDIDAS

if mod(kt,Tmed) == 0 && kt ~= 0 %Sólo incorporamos una medida cada 3 horas

    z = variables.zmed(med); %Tomamos la medida (observaciones)

    r = sqrt(xp(1)^2 + xp(2)^2); %Recalculamos el Radio
    
    %Creamos el vector de transformación de medidas a partir del modelo de
    %cómo se generan las medidas
    H = [xp(1)/r xp(2)/r 0 0];
    
    %Redimensionamos la matriz P (forma matricial 4x4)
    P = reshape(xp(5:20),[4, 4]);

    %Ganancia de Kalman
    K_gain = P * H.' * inv(H * P * H.' + R) ;
    
    % Estimación corregida por la ganancia de Kalman y un factor que
    % contiene la diferencia entre medidas y medidas estimadas
    xp(1:4) = transpose(xp(1:4)) + K_gain*(z - r);
    
    % Se corrige también la estimación del error
    P = P - (K_gain * H * P);
    
    p_array = reshape(P, [1, 16]);

    xp = [xp(1) xp(2) xp(3) xp(4) p_array];
    
    % Radio observado
    r_est_medidas(k) = z; 
    
    med = med + 1; %Incrementamos el iterador (siguiente medida)

end

% EVALUACIÓN DE RESULTADOS

%ESTO AQUI HAY QUE GUARDARSE MAS COSAS OSEA LAS 20 PQ TE HACE FALTA PARA EL
%ERROR. GUARDATE MAS COSAS DE XP Y MIRA EL GUIÓN.

%Almacenamos la distancia estimada para cada iteración
r_est(k) = sqrt(xp(1)^2 + xp(2)^2);

%Realizamos la comparación con las medidas reales y lo almacenamos.
%Error cometido en la posición
error_2d = variables.Cxy_true(k,:) - xp(1:2);
error_cometido(k) = sqrt( error_2d(1)^2 + error_2d(2)^2 ); %Expresamos el error en formato cuadrático

% Error cometido estimado
sigmax = P(1,1); % Elemento 1 de la diagonal 
sigmay = P(2,2); % Elemento 2 de la diagonal
error_estimado(k) = sqrt(sigmax^2 + sigmay^2);

orbita(:,k) = xp(1:2); %Guardamos el valor estimado de la órbita

r_orig(k) = sqrt(variables.Cxy_true(k, 1)^2 + variables.Cxy_true(k,2)^2); 
%Almacenamos el valor verdadero de la órbita para comparación


k = k + 1;
end

%Comparar la desviación estándar del ruido de medida con la desviación
%estándar del error entre rk y rk_est
% Radio original

r_est_medidas(r_est_medidas == 0) = NaN;
figure(1);
plot(r_est_medidas,'*');
hold on;
plot(r_est);
hold on;
plot(r_orig)
title('Radio medido');
legend('medidas de radio','radio estimado','radio verdadero')

% Comparación de las órbitas
figure(2);
plot(orbita(1,:),orbita(2,:))
hold on
plot(variables.Cxy_true(:, 1), variables.Cxy_true(:, 2))
title('Comparación de órbitas')
legend('Órbita estimada', 'Órbita real');

% Comparación del error cometido y el error estimado
figure(3);
plot(error_cometido);
hold on
plot(error_estimado);
ylabel('Error');
xlabel('Instante de medida');
legend('Error cometido', 'Error estimado');
title('Comparación errores cometido y estimado');




