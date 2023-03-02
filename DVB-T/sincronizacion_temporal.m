
function  muestra_inicial=sincronizacion_temporal(x,NFFT,Ng)

% Esta función determina cual es la muestra inicial en la que comienza el
% primer símbolo de una señal OFDM. Para ello hace uso del número de
% portadoras teóricas y de la duración del intervalo de guarda en muestras
%
% x: Señal OFDM
% NFFT: Número de portadoras teóricas (debe ser 2048 u 8192)
% Ng: Duración del intervalo de guarda (en muestras)
%
% muestra_inicial: Muestra en la que comienza el primer símbolo OFDM

Nd=4*NFFT;          % Número suficiente de símbolos OFDM

% En este caso hacemos uso de la funcion maximum correlation metric pero
% usando como Nw la duración del intervalo de guarda.

for n=1:Nd
    
    suma=0;
    for i=1:Ng
        suma=suma+conj(x(n+i))*x(n+i+NFFT);
    end
    p(n)=suma/Ng;
end

% El primer máximo de la funcion p(n) coincide con la muestra en la que
% comienza el primer símbolo OFDM. Dado que la función p(n) presenta mucho
% rizado cerca de los picos máximos nos quedamos unicamente con las 2048
% primeras muestras (si la señal es 2K) o con las 8192 primeras (si la
% señal es 8K) ya que en ellas se encontrará el primer máximo y no nos
% confundiremos con otro posterior.

[~,muestra_inicial]=max(p(1:NFFT));

% Se representa la funcion maximum correlation metric 

figure(8)
plot(abs(p))
title('Sincronizacion temporal')
ylabel('|pNFFT,Nw(n)|')
xlabel('n')
grid on

end

