function [d] = static_filter(coef, z)
% Función que realiza el filtrado fijo. Se
% produce un filtrado con los coeficientes dados mediante la función
% "circshift()" con la cual retrasamos una muestra la señal z para poder ir
% filtrando muestra a muestra.
% ENTRADA:
% - coef : Vector que contiene los coeficientes del filtro
% - z: señal que se va a filtrar
% SALIDA:
% - d : señal filtrada

M = length(coef);
window = zeros(1, M);
d = zeros(1, length(z));

for n=1:length(z)
    window = circshift(window, 1); % Se va desplazando la ventana
    window(1) = z(n); 
    d(n) = window * coef; %Realizamos el filtrado
end

