function [W] = pesos(tn)
%PESOS: Función que calcula los pesos del beamformer para cada f

N = length(tn); %longitud del vector de retardos
flim = 129;     %Barrido de frecuencias
W = zeros(flim,N); %Inicializamos el vector de pesos

    for f = 1:flim       
        for i = 1:N

            W(f,i) = (1/N)*exp(j*2*pi*tn(i)*f); %Vector 1x7

        end
    end
end