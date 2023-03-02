%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% DVB-T RECEIVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autores: Javier Lobato Martín y Sergio Zapata Caparrós

%% Visualiza la PSD de una señal TDT
clc
clear all
close all
%fichero='tdt_486.int16'; Fc=486e6; Fs=20e6;
fichero='tdt_482.int16'; Fc=482e6; Fs=10e6;
%fichero='tdt_490.int16'; Fc=490e6; Fs=10e6;
%fichero='tdt_506.int16'; Fc=506e6; Fs=10e6;
%fichero='tdt_514.int16'; Fc=514e6; Fs=10e6;
%fichero='tdt_634.int16'; Fc=634e6; Fs=20e6;
%fichero='tdt_650.int16'; Fc=650e6; Fs=10e6;

% Carga de los datos
xx=load_bbfile(fichero,'int16');
Nfft=2048;
x=xx(1:100000); 


% %Filtro paso baja
% [A, B] = butter(4, 4e6/(Fs/2)); %[b,a] = butter(nth,fc/(fs/2));
% % freqz(A, B, [], Fs); %Comprobación del filtro
% x = filter(B, A, x);

figure(1);
[Pxx, Fxx] = pwelch(x,4096,2048,4096,Fs, 'centered','power');
plot(Fxx, 10*log10(Pxx));
title('Canal');
xlim([-Fs/2 Fs/2]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');

% subplot(2,1,2);
% [S,F,T]=specgram_complex(x,Nfft,Fs,hamming(Nfft));
% S=10*log10(abs(S));
% Fv=(F+Fc)/1e6;
% Tv=1000*T;
% subplot(2,1,2);
% imagesc(Fv,Tv,S'); colormap(jet);
% xlabel('Freq(MHz)');
% ylabel('Time(ms)');


%% Parámetros recepción
Ts = (7/64)*1e-3; % Periodo de símbolo
T = 1/8e6; %periodo de muestreo 
Fs2 = 8e6;     % 8 MHz de nuestra segunda frecuencia de muestreo
Nss = 7168;    % Número de muestras por símbolo
Tsu = 896e-6; % Tiempo de símbolo útil



x1 = resample(x, 4,5); % Remuestreamos a 8 MHz
%Ahora tenemos 7168 muestras/símbolo

%%

%Aquí podemos comprobar todas las portadoras que tenemos, por cada pico
%vamos a tener una portadora

x1freq = fft(x1,8192); %Haciendo la FFT decodificamos la señal con las 8192 portadoras separadas delta f
x1freq = fftshift(x1freq);
figure(2)
plot(abs(x1freq))
title('Disposición de las portadoras de la señal')

%% BARRIDO PARA OBTENER EL TAMAÑO DE LA BANDA DE GUARDA


tam=[28, 56, 112, 224]*1e-6*Fs2;% Número de muestras de intervalo de guarda
%1/32 1/16 1/8 1/4

correla = zeros(4,(Ts*Fs2)+tam(4)); %4 posibles intervalos de guarda y la longitud en muestras de 1 simbolo
%El resto de valores se quedarán en 0, ya que para los tres primeros casos
%no se evalúa hasta esa longitud

l_symb = Nss; %Nss es el número de muestras que vamos a tener en un símbolo

for i=1:4
    
    i_guarda = tam(i);

    for j = 1:(l_symb+tam(i)) %Número de muestras que corresponden a un símbolo + intervalo de guarda correspondiente
    
       correla(i,j) = max(abs(xcorr(x1(j:j+i_guarda),x1(j+l_symb:j+i_guarda+l_symb),0)));         

    end

end
%EXPLICACIÓN: Para cada valor de intervalo de guarda, evaluamos la
%correlación que hay entre una porción de la señal y la misma porción un
%símbolo despues. Esto se evalúa un número de veces igual a la longitud de
%un símbolo, ya que vamos probando para cada muestra del total de un símbolo.
%Al ejecutarse el bucle 4 veces, obtendremos 4 líneas con el valor de la 
%correlación para cada elemento. 




%Hacemos un plot que nos muestre todas las funciones de correlación.
%Podemos comprobar como en la última tenemos la función de correlación
%esperada
figure(3);
plot(correla(1,:));
figure(4);
plot(correla(2,:));
figure(5);
plot(correla(3,:));
figure(6);
plot(correla(4,:));


%A continuación, usando max otra vez sacamos el máximo de correla, lo que
%nos va a dar el tamaño del intervalo de guarda y el desplazamiento justo
%que nos da el comienzo de del símbolo.
maxfila = max(correla,[],2);
[~,fila_max] = max(maxfila);
i_guarda = tam(fila_max);


[maximo,ind_max] = max(correla(fila_max,:));

muestra_primer_simbolo = sincronizacion_temporal(x1, Nss, i_guarda); %Estos parámetros están adaptados para 8 MHz de frecuencia de muestreo.
%Corresponden al modo 8K y a la duración del intervalo de guarda (1/4)
%%
x1_syncro = x1(muestra_primer_simbolo:end);

%A continuación, hacemos un reshape para reorganizar las muestras en forma
%de matriz, para que adquieran la disposición de un símbolo por columna. De
%esta manera, eliminaremos los intervalos de guarda y practicaremos la FFT.

len_simbolo = Nss+i_guarda;
n_simb = floor(length(x1_syncro)/len_simbolo);

%Le quitamos a x1 el número de muestras por el final que sobran, para que
%reshape funcione correctamente
resto = mod(length(x1_syncro),len_simbolo);
x1_syncro = x1_syncro(1:n_simb*len_simbolo);

%Hacemos reshape
X_MAT = reshape(x1_syncro, [len_simbolo, n_simb]);

%Eliminamos los intervalos temporales de guarda, que corresponde con
%eliminar las primeras 1792 filas de nuestra matriz:
X_MAT(1:i_guarda,:) = [];


%A continuación, practicamos la fft, que demodula en una sola operación
%todas las portadoras a la vez. Usaremos la dimensión = 7168 para
%adaptarnos al tamaño de nuestras portadoras.

X_MAT_DEMOD = fft(X_MAT,Nss); %Nss es 7168
X_MAT_DEMOD_SHIFT = fftshift(X_MAT_DEMOD);

%Para comprobar que estamos realizando la demodulación de manera correcta,
%promediamos los valores de todos los símbolos para cada fila. Una fila
%corresponde a una portadora. Lo representamos para comprobar las
%portadoras y sus potencias correspondientes

X_MAT_DEMOD_MEAN = mean(abs(X_MAT_DEMOD_SHIFT),2);

figure(7)
plot(X_MAT_DEMOD_MEAN);
title('Disposición de las portadoras en la señal demodulada')
%Como podemos observar, tenemos la señal demodulada correctamente. Podemos
%comprobar fácilmente que la demodulación es correcta comprobando las
%portadoras continuas (aquellas que se envían con una potencia mayor). La
%primera y la última son portadoras continuas, por ejemplo. 

%% DEMODULACIÓN DE LAS PORTADORAS TPS

%A continuación, observamos en el estándar dónde se encuentran las
%portadoras TPS y realizamos una demodulación DBSPK.
%Debido a que estamos en el modo 8K, según el estándar hay 68 portadoras
%TPS repartidas en las portadoras, las primeras se encuentran en las
%portadoras: 13, 50, 209, 346, etc
%Para comprobar la demodulación bastará con mostrar que la nube de puntos
%corresponde con una de DBPSK. 

%Primero tomamos las portadoras útiles, que comienzan en la 178 y terminan
%en la 6994
PORTADORAS_UTILES = X_MAT_DEMOD(178:6994,:);

TPS = PORTADORAS_UTILES(34,:);

%Para realizar la demodulación DBSPK tenemos que obtener la diferencia de
%fase que hay entre símbolos. Podemos realizar esto multiplicando cada
%símbolo por el conjugado del anterior. Utilizamos el conjugado porque
%queremos ver la diferencia de fase, con el conjugado obtenemos
%la fase de signo opuesto para cada símbolo.

TPS_conj = conj(TPS);

TPS_demod = TPS(2:end).*(TPS_conj(1:end-1));

figure(9)
plot(real(TPS_demod),imag(TPS_demod),'b.')
title('Nube de puntos TPS');
grid("on")

