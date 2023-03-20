function [d] = adaptable_filter(coef, z, mu)

% Función que realiza el filtrado adaptable. Para la primera iteración se
% produce un filtrado con los coeficientes iniciales mediante la función
% "circshift()" con la cual retrasamos una muestra la señal z para poder ir
% filtrando muestra a muestra.
% A partir de la segunda iteración, se calcula el error cometido y se
% aplica descenso en gradiente para adaptar los coeficientes del filtro.

% ENTRADA:
% - coef : Vector que contiene los coeficientes del filtro
% - z: señal que se va a filtrar
% - mu: Paso descenso en gradiente
% SALIDA:
% - d : señal filtrada

%Inicializamos variables
M = length(coef);
window = zeros(1, M);
d = zeros(1, length(z));


for n=1:length(z)
    
    window = circshift(window, 1); % Se va desplazando la ventana
    window(1) = z(n); 
    
    d(n) = window * coef; %Realizamos el filtrado

    simb = qamdemod(d(n),16); %Obtenemos el símbolo al que corresponde s(n-u)
    simb_mod = qammod(simb,16); %Esto nos da s(n)
    e = simb_mod - d(n); %Calculamos el error cometido
    z_mod = (conj(window) * transpose(window)); %Obtenemos la estimación de la energía de la señal
    alpha = (mu/z_mod)*e;
    coef = coef + alpha.*(window'); %Actualizamos los nuevos coeficientes del filtro

end

