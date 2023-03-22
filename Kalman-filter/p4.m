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

w = [0 0 sig2w sig2w];   % Vector de las perturbaciones
Q = diag(w);                % Matriz Covarianza(de las perturbaciones)

P = diag([var_pos var_pos var_v var_v]); % Matriz de errores (diagonal, inicial)
%R = diag([0.1 0.1]); % Matriz de errores de medida (diagonal)
R = 0.1; % Error de medida

%% Construcción del filtro

%Primero de todo, conformamos el vector xp inicial, concatenación de x y los
%elementos de la matriz P reorganizados:
xp = [x reshape(P,[1,16])];
orbita = zeros(2,length(variables.tt));
med = 1; %Contador que nos va a indicar qué medida tenemos que tomar
k = 1;   %Contador que nos va a indicar en qué instante nos encontramos (en formato entero)
error = zeros(1,length(variables.tt)); %Almacenaremos el error en las 240 iteraciones
pos = zeros(2,length(variables.tt));   %Almacenaremos la posición estimada en las 240 iteraciones
r_est = zeros(1,length(variables.tt)); %Almacenaremos el radio estimado de órbita
r_est_medidas = zeros(1,240); %Aquí almacenaremos las medidas tomadas cada 3 horas (para comparar luego)
r_orig = zeros(1, 240); %Aquí vamos a almacenar el radio verdadero


for kt = 0:Tpred:24-Tpred


% PREDICCIÓN
%Para la primera iteración, tenemos ya conformado el vector xp incial. Para
%el resto de iteraciones, tenemos el vector xp de salida de la iteración
%anterior

%r = sqrt(xp(1)^2 + xp(2)^2); % Radio

%CORRIGE LA MATRIZ DE TRANSICIÓN DE LA FUNCIÓN DIFEQ
%Debemos conformar la matriz de transición F
% elemento31 = GM * (3*xp(1)^2 - r^2) / r^5;
% elemento32 = GM * (3*xp(1) * xp(2)) / r^5;
% elemento41 = GM * (3*xp(1) * xp(2)) / r^5;
% elemento42 = GM * (3*xp(2)^2 - r^2) / r^5;
% F = [0 0 1 0; 0 0 0 1; elemento31 elemento32 0 0; elemento41 elemento42 0 0];

%OJO ESTO NO SE SI VA AQUI (creo que si): Obtenemos la estima del error
%P = F * P * F' + Q;
%P = F * P + P * F' + Q;


%COMPROBAR SI ESTO ES ASÍ
tini = kt;
tfin = tini + Tpred;

%Resolvemos con ode45 las ecuaciones diferenciales
[~, xp] = ode45(@difeq, [tini, tfin], xp);


xp = xp(end, :); % Tomamos las 4 variables de estado


% INCORPORACIÓN DE MEDIDAS

%COMPROBAR TODO ESTO, NO SE SI ESTÁ BIEN
if mod(kt,Tmed) == 0 && kt ~= 0 %Sólo incorporamos una medida cada 3 horas

    z = variables.zmed(med); %Tomamos la medida (observaciones)

    r = sqrt(xp(1)^2 + xp(2)^2); %Recalculamos el Radio
    
    H = [xp(1)/r xp(2)/r 0 0]; %Creamos el vector de transformación de medidas
    
    P = reshape(xp(5:20),[4, 4]);

    %Ganancia de Kalman
    K_gain = P * H.' * inv(H * P * H.' + R) ; %Obtenemos la ganancia de Kalman

    %xp(1:4) = transpose(xp(1:4)) + K_gain*(z - H*transpose(xp(1:4))); %Corregimos nuestra estimación con las medidas obtenidas
    xp(1:4) = transpose(xp(1:4)) + K_gain*(z - r);

    P = P - (K_gain * H * P); %Corregimos nuestra estimación del error con las medidas obtenidas
    

    dp_array = reshape(P, [1, 16]);

    % dxp tiene que ser un vector columna(formato de la función ode45)
    xp = [xp(1) xp(2) xp(3) xp(4) dp_array];

    r_est_medidas(k) = z; 
    
    med = med + 1; %Incrementamos el iterador

end

% EVALUACIÓN DE RESULTADOS

%ESTO AQUI HAY QUE GUARDARSE MAS COSAS OSEA LAS 20 PQ TE HACE FALTA PARA EL
%ERROR. GUARDATE MAS COSAS DE XP Y MIRA EL GUIÓN.

%Almacenamos la distancia estimada para cada iteración
r_est(k) = sqrt(xp(1)^2 + xp(2)^2);

%Realizamos la comparación con las medidas reales y lo almacenamos.
error_2d = variables.Cxy_true(k,:) - xp(1:2);
error(k) = sqrt( error_2d(1)^2 + error_2d(2)^2 ); %Expresamos el error en formato cuadrático

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

% HAY QUE TRAZAR LA ÓRBITA ESTIMADA DEL SATÉLITE (con el radio)
figure(2);
plot(orbita(1,:),orbita(2,:))
hold on
plot(variables.Cxy_true(:, 1), variables.Cxy_true(:, 2))
title('Comparación de órbitas')



