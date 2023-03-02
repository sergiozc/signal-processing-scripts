
function  muestra_inicial=sincronizacion_temporal(x,NFFT,Ng)

% Esta funci�n determina cual es la muestra inicial en la que comienza el
% primer s�mbolo de una se�al OFDM. Para ello hace uso del n�mero de
% portadoras te�ricas y de la duraci�n del intervalo de guarda en muestras
%
% x: Se�al OFDM
% NFFT: N�mero de portadoras te�ricas (debe ser 2048 u 8192)
% Ng: Duraci�n del intervalo de guarda (en muestras)
%
% muestra_inicial: Muestra en la que comienza el primer s�mbolo OFDM

Nd=4*NFFT;          % N�mero suficiente de s�mbolos OFDM

% En este caso hacemos uso de la funcion maximum correlation metric pero
% usando como Nw la duraci�n del intervalo de guarda.

for n=1:Nd
    
    suma=0;
    for i=1:Ng
        suma=suma+conj(x(n+i))*x(n+i+NFFT);
    end
    p(n)=suma/Ng;
end

% El primer m�ximo de la funcion p(n) coincide con la muestra en la que
% comienza el primer s�mbolo OFDM. Dado que la funci�n p(n) presenta mucho
% rizado cerca de los picos m�ximos nos quedamos unicamente con las 2048
% primeras muestras (si la se�al es 2K) o con las 8192 primeras (si la
% se�al es 8K) ya que en ellas se encontrar� el primer m�ximo y no nos
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

