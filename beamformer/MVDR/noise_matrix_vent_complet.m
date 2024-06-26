function [corr_noise] = noise_matrix_vent_complet(N, freq, win, Ltrama, Lfft, muestras_ruido, y)

noise = y(1:muestras_ruido, :);

% Garantizamos que el número de muestras del ruido sea divisible en tramas 
[m,~]=size(noise);
resto=mod(m,Ltrama/2);
noise=noise(1:m-resto,:);

% Se obtiene el número de muestras que tendrá el ruido sobre el que se
% cálcula la matriz de correlación espacial
[m,~]=size(noise);
Ntramas=2*(m/Ltrama)-1;


% Matriz de NxN para cada frecuencia
corr_noise=zeros(N,N,Lfft/2 +1);
iter = 1;
for ntrama=1:Ntramas
    trama_noise = noise(iter:iter + Ltrama ,:); %Tomamos la porción de señal del canal correspondiente
    trama_f = fft(win.*trama_noise, Lfft);
    trama_f = trama_f(1:Lfft/2+1,:);

    for i=1:N % Sensor i
        for j=1:N % Sensor j
            for k=1:length(freq) % Frecuencia k
                corr_noise(i,j,k) = corr_noise(i,j,k) + trama_f(k,i) * trama_f(k,j)';
            end
        end
    end
    iter = iter + 127; %Actualizamos el iterador
end

corr_noise = corr_noise / (Ntramas); % Normalización
end

