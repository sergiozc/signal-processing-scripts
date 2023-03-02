
"""
ANTONIO SIMÓN MARTÍN
SERGIO ZAPATA CAPARRÓS

"""

import numpy as np
import matplotlib.pyplot as plt
from skimage import io, color, util
import scipy.ndimage as ndi
from scipy import signal


def TF_imagen(I):
    
    '''Función que realiza la transformada de fourier en dos dimensiones
    dada una imagen'''
    
    # Transformada de Fourier y espectro
    # Cargar imagen y convertir a niveles de gris
    I = io.imread('mandrill.tif',as_gray=True)
    plt.figure()
    plt.subplot(1,2,1), plt.imshow(I,cmap='gray')
    
    
    # Se selecciona una subimagen con una textura clara
    B = I[130:130+100,150:150+200] #Sección seleccionada
    plt.subplot(1,2,2), plt.imshow(B,cmap='gray')
    N,M = B.shape
    imsize = N*M
    plt.figure()
    plt.imshow(B,cmap='gray')
    
    
    # Computar FFT2
    
    F = np.fft.fft2(B) #Realizamos la transformada en las filas y la transformada en las columnas
    
    # Usar escala logaritmica para correcta visualizacion
    plt.imshow(10*np.log10(np.abs(F)**2/imsize),cmap='hot')
    plt.colorbar()
    plt.title('Transformada de Fourier de la imagen')
    
    # Cuadrantes reordenados (frec. (0,0) en el centro)
    
    F=np.fft.fftshift(F) #Desplazamos para que esté centrado
    plt.figure()
    plt.imshow(10*np.log10(np.abs(F)**2/imsize),cmap='hot')
    plt.colorbar()
    plt.title('Transformada de Fourier desplazada de la imagen')
    
    pass

def HEQ(I):
    
    '''Función que realiza la comparación tanto de histogramas como de imágenes
    mediante una reorganización de los componentes de intensidad, dada una imagen'''
    
    plt.figure() # Mostrar imagen original
    plt.imshow(I,cmap = 'gray')
    plt.title('Imagen original')
    plt.axis('off')
    Nbins = 256 #Numero de subintervalos
    bins=np.arange(-0.5,Nbins) # Definir limites de los subintervalos
    Hist,bins=np.histogram(I,bins) # Cómputo del histograma
    centers = 0.5*(bins[1:]+bins[:-1]) # Calcular centros de los histogramas
    plt.figure() # Dibujar histograma original
    plt.stem(centers,Hist)
    plt.title('Histograma Original')

    #Este histograma muestra la distribución de intensidades dentro de la imagen original
    #observamos que las componentes de intensidad están todas muy juntas entre ellas 
    
    T = Hist.cumsum() / Hist.sum() * 255
    T.round().astype(int)
    I = T[I]
    plt.figure() #Imagen realzada
    plt.imshow(I,cmap='gray')
    plt.title('Imagen realzada')
    plt.axis('off')
    
    
    Hist2,bins=np.histogram(I,bins) # Cómputo del histograma
    plt.figure() # Dibujar histograma original
    plt.stem(centers,Hist2)
    plt.title('Histograma Resultante')
    
    pass

def filtrado_lineal(I,J):
    
    '''Funcion que implementa el filtrado lineal dadas la imagen original y
    la imagen degradada y devuelve la imagen restaurada
    
    ENTRADA: 
        I: Imagen original
        J: Imagen degradada
        
    SALIDA:
        R: Imagen restaurada
        h: Respuesta impulsiva
        H: Respuesta en frecuencia
        
    '''
    
    #FILTRADO LINEAL
    #Ahora hay que operar con el filtro de 1 dimensión
    # Generar respuesta en frecuencia deseada Hd
    Lfir=11 # Orden FIR 1D (impar)
    fc=0.18 # Frecuencia de corte 1D (f=1 --> frecuencia de muestreo)
    Lcut=int(np.round(fc*Lfir)) # bin de frecuencia para fc
    lcent=int((Lfir-1)/2) # punto central (freq. cero)
    Hd_1D=np.zeros(Lfir) # filtro Horizontal
    Hd_1D[lcent-Lcut:lcent+Lcut+1]=np.ones(2*Lcut+1)
    
    H1d_shift = np.fft.ifftshift(Hd_1D) #Reorganizamos los cuadrantes de frecuencia
    
    #Ahora vamos a introducir una fase lineal para que quede centrado en lcent
    w=2*np.pi*np.linspace(0,(Lfir-1),Lfir)/Lfir
    H1d_lineal = np.abs(H1d_shift) * np.exp(-1j * w * lcent)
    
    Hd = np.outer(H1d_lineal,H1d_lineal) # filtro 2D separable cuadrado (filtros H y V identicos)
    
    h = np.real(np.fft.ifft2(Hd)) #he cambiado real por abs
    
    H = np.fft.fft2(h,(1024,1024))
    H = np.fft.fftshift(H)
    #r_freq_real = np.real(r_freq)
    
    #Ahora realizamos el filtrado de la imagen contaminada con el filtro
    R = ndi.convolve(J,h)
    
    return R, h, H

def filtrado_no_lineal(I, J, N):
    
    '''Funcion que implementa el filtrado no lineal
    mediante el filtro de mediana
    
    ENTRADA: 
        I: Imagen original
        J: Imagen degradada
        N: Tamaño del filtro de mediana
    SALIDA:
        R_median: Imagen restaurada
    '''

    R_median = ndi.median_filter(J, size = 3)
    
    return R_median

'''    
#Cargamos y visualizamos la imagen en formato Bitmap (.bmo)
I=io.imread('lena.bmp') #Contiene 3 imagenes con componentes R, G y B
plt.figure()
plt.imshow(I)
plt.axis('off')

#Podemos convertir a escala de grises
I_=io.imread('lena.bmp',as_gray=True)#Ahora I contiene una solo imagen, una matriz
plt.figure()
plt.imshow(I_, cmap='gray')
plt.axis('off')
'''

#%% TECNICA HEQ: trata de redistribuir las intensidades de los pixels de forma que el histograma cmbia a otro 
# predeterminado, la ecualizacion de histogramas (HEQ) tiene como objetivo aumentar el contraste de la imagen

#Histograma: representa la distribucion de la intendad luminosa de una imagen. Se obtiene dividiento el intervalo
#de intensidades en un determinado numero de subintervalos Nbins
#Para cada subintervalo se representa el numero de pixels cuya intensidad ce dentro del mismo

I=io.imread('tire.tif') # Lectura de imagen tiff

HEQ(I)

#Apreciamos cómo se realza la imagen resultante respecto a la original y podemos ver
#también que las componentes de las intensidades se han redistribuido en el histograma
#resultante, en comparación con el histograma original


#%% TRANSFORMADA DE FOURIER 2D: Realizamos la transformada 1D a las filas de la imagen 
#y seguidamente a las columnas

I = io.imread('mandrill.tif',as_gray=True)

TF_imagen(I)

#Observamos como, principalmente, tenemos componentes en baja frecuencia, ya que
#el espectro se basa en componentes de 20 en la escala de frecuencia mostrada. También
#podemos apreciar la eficiencia de desplazar la transformada ya que hemos conseguido centrar
#los componentes de mayor frecuencia en el centro

#%%EJEMPLO FILTRADO

# Generar respuesta en frecuencia deseada Hd
Lfir=11 # Orden FIR 1D (impar)
fc=0.18 # Frecuencia de corte 1D (f=1 --> frecuencia de muestreo)
Lcut=int(np.round(fc*Lfir)) # bin de frecuencia para fc
lcent=int((Lfir-1)/2) # punto central (freq. cero)
Hd_1D=np.zeros(Lfir) # filtro Horizontal
Hd_1D[lcent-Lcut:lcent+Lcut+1]=np.ones(2*Lcut+1)
Hd=np.outer(Hd_1D,Hd_1D) # filtro 2D separable cuadrado (filtros H y V identicos)
plt.figure()
plt.imshow(Hd)
plt.colorbar()
plt.title('Respuesta en frecuencia deseada')

#%% RESTAURACION DE IMAGENES: El objetivo es que a partir de una imagen distorsionada
#obtiene otra imagen que represente lo mejor posible la imagen original


#FILTRADO LINEAL:
I=io.imread('lena.bmp',as_gray=True) # Cargar imagen
J=util.random_noise(I, mode='s&p') #Contaminada con ruido sal y pimienta

#Realizamos un filtrado lineal paso-baja
R, h, r_freq=filtrado_lineal(I,J)


#Representamos las tres imagenes
figure, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
ax1.imshow(I, cmap='gray')
ax1.set_axis_off()
ax1.set_title('Imagen original')
ax2.imshow(J, cmap='gray')
ax2.set_axis_off()
ax2.set_title('Imagen contaminada')
ax3.imshow(R, cmap='gray')
ax3.set_axis_off()
ax3.set_title('Imagen restaurada (filtro lineal)')
plt.tight_layout()

#Apreciamos que los contornos de la imagen no se recuperan muy bien,
#debido a la característica paso-baja y no se consigue totalmente eliminar
#el ruido sal y pimienta, aun así se ha mejorado tras el filtrado

plt.figure(2)
plt.plot(h)
plt.title('Respuesta al impulso') 

plt.figure(3)
plt.imshow(np.abs(r_freq))
plt.colorbar()
plt.title('Respuesta en frecuencia')

#Comprobamos que se corresponde a un filtro paso-baja 2D  y que su respuesta
#en frecuencia es muy similar a la deseada

#%%FILTRADO NO LINEAL:

I=io.imread('lena.bmp',as_gray=True) # Cargar imagen
J=util.random_noise(I, mode='s&p') #Contaminada con ruido sal y pimienta

N=3
#N=5
#N=7

R_median = filtrado_no_lineal(I, J, N)


#Representamos las tres imagenes
figure, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
ax1.imshow(I, cmap='gray')
ax1.set_axis_off()
ax1.set_title('Imagen original')
ax2.imshow(J, cmap='gray')
ax2.set_axis_off()
ax2.set_title('Imagen contaminada')
ax3.imshow(R_median, cmap='gray')
ax3.set_axis_off()
ax3.set_title('Imagen restaurada filtro mediana 3x3')
plt.tight_layout()

#Con un tamaño de radio de 3x3 se restaura la señal correctamente con los
#bordes bastante nítidos
#Observamos que al aumentar el entorno de radio donde trabaja la mediana
#la imagen restaurada es peor, ya que se traduce en realizar la mediana
#de más número distintos





