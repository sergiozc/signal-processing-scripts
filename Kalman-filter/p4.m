%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PRÁCTICA 4: Filtro de Kalman %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autores: Javier Lobato Martín y Sergio Zapata Caparrós
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all; 
close all;

dist = 35786;   % Distancia en km a la superficie terrestre
v = 11068;      % Velocidad orbital (km/h)
G = 6.6742e-11; % Cte gravitacional (m3/kg*s2)
M = 5.9736e24;  % Masa de la Tierra (kg)

x = [42164 0 0 11068];         % Vector de estado [x y x(.) y(.)]
Tpred = 0.1;    % Tiempo de predicción (0.1 h)
Tmed = 3;       % Tiempo de medidas (3 h)
var_pos = 50;   % Varianza de error posición (km2)
var_v = 1;      % Varianza de error velocidad (km2/h2)
global Q        % Matriz covarianza error
global GM;



variables = load('Variables_Orbita.mat');
sig2w = variables.sig2w;    % Perturbación de la aceleración
GM = variables.GM;

w = [0; 0; sig2w; sig2w];   % Matriz de las perturbaciones
Q = diag(w);                % Matriz Covarianza

for t = 0:Tpred:24-Tpred

P = diag([var_pos var_pos var_v var_v]); % Matriz de errores diagonal

% ... PREDICCIÓN
% Hay que calcular xp que es la concatenacion de x y P
[tp, xp] = ode45(@difeq, [tini, tfin], xp);

% Estas son las 4 variables de estado
xp = xp(end, :); %Hay que comparar estos valores con los reales

% INCORPORACIÓN DE MEDIDAS
% ¿Hay medidas?


end
