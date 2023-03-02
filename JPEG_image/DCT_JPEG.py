'''Implementación de un codificador JPEG de imagen con transformada DCT 2D'''

'''Autores: Antonio Simón Martín y Sergio Zapata Caparrós'''


import numpy as np
import matplotlib.pyplot as plt
from skimage import io, color
from scipy.fft import dctn, idctn

def psnr(I, I_reconstruida):
    
    N1=len(I)
    N2=len(I)
    
    peak = np.max(np.abs(I))
    denom = (1/float(N1*N2))*np.sum((I-I_reconstruida)**2)
    psnr=10.0*np.log10((peak**2)/denom)
    
    return psnr

def busca_ceros(B):
    
    ceros = 0
    
    for i in range (8):
        for j in range (8):
            if B[i][j] == 0:
                ceros = ceros + 1
            else:
                pass
            
    return ceros



plt.close('all')



#En primer lugar leemos la imagen y la convertimos en escala de grises
#Le restamos la media para eliminar el componente en continua (rango [-127,128])
I = io.imread('lena.bmp', as_gray=True).astype(np.double)*255 -128

#Si las componentes de la imagen no se pueden dividir bloques de 8, añadimos ceros
resto = len(I) % 8 
if resto == 0:
    pass
else:
    nceros = 8 - resto
    ceros = np.zeros(nceros)
    I = np.append(I,ceros, axis = 0)
    I = np.append(I, ceros, axis = 1)


mask = np.zeros((8, 8))
I_reconstruida = np.zeros_like(I)
bloque = np.zeros((8,8)) #Combinación lineal de las 64 frecuencias

for i in range (64):
    for j in range(64):
        #Nos quedamos con bloques 8x8
        bloque = I[i*8:(i*8)+8, j*8:(j*8)+8]
        
        #Hacemos la DCT a cada valor del bloque
        bloque_dct = dctn(bloque, norm='ortho')
        
        #Ahora tenemos que quedarnos con las frecuencias mas bajas, es decir 
        #con un bloque 4x4 de los coeficientes calculados
        
        N = 4 #Factor de muestreo zonal
        
        mask[:N, :N] = 1 #Creamos una matriz de zeros y ponemos a 1 las componentes
                         #con las que nos queremos quedar.
                       
        bloque_muestreo = bloque_dct * mask #Ponemos a cero las componentes que no nos interesan
         
        #Reconstruimos la señal
        I_reconstruida[i*8:(i*8)+8, j*8:(j*8)+8] = idctn(bloque_muestreo, norm='ortho') #Haciendo el muestreo zonal
        

#Hacemos la DCT 2D
D = dctn(I, norm='ortho')


figure, (ax1, ax2)=plt.subplots(1, 2, figsize=(10, 5))
ax1.imshow(I, cmap='gray')
ax1.set_title('Imagen original')
ax2.imshow(np.log(D**2))
ax2.set_title('DCT 2D')
plt.tight_layout()
#Apreciamos en la esquina superior izquierda de la imagen de la transformada,
#la componente en continua, y conforme nos alejamos, componentes con mayor frecuencia

figure, (ax1, ax2, ax3)=plt.subplots(1, 3, figsize=(15, 5))
ax1.imshow(I, cmap='gray')
ax2.set_axis_off()
ax1.set_title('Imagen original')
ax2.imshow(I_reconstruida, cmap='gray')
ax2.set_axis_off()
ax2.set_title('Imagen reconstruida')
ax3.imshow(np.abs(I-I_reconstruida), cmap='gray')
ax2.set_axis_off()
ax3.set_title('Error')
plt.tight_layout()
#Observamos que, al quedarnos con los primeros coeficientes (coeficientes de bajas frecuencias),
#el error se produce en los contornos de la imagen, es decir, en las componentes con
#mayor frecuencia. Cambiando el factor de muestreo zonal, comprobamos como varía el error.

#Calculamos ahora la peak signal to noise ratio, medida para la calidad de la señal empleada en imágenes

PSNR = psnr(I,I_reconstruida)

print('La PNSR es de',PSNR, 'dB')


#%% CODIFICADOR JPEG: Con cuantizador

#En primer lugar leemos la imagen y la convertimos en escala de grises
I_JPEG = io.imread('lena.bmp', as_gray=True).astype(np.double)*255-128

#Si las componentes de la imagen no se pueden dividir bloques de 8, añadimos ceros
resto = len(I_JPEG) % 8 
if resto == 0:
    pass
else:
    nceros = 8 - resto
    ceros = np.zeros(nceros)
    I_JPEG = np.append(I_JPEG,ceros, axis = 0)
    I_JPEG = np.append(I_JPEG, ceros, axis = 1)



#Inicializamos variables
q0=np.zeros(10)
ind=0
PSNR_JPEG=np.zeros(10)
I_JPEG_reconstruida = np.zeros_like(I_JPEG)

for q in range (5, 55, 5):
    
    q0[ind] = q
    bloque = np.zeros((8,8)) #Combinación lineal de las 64 frecuencias
    
    
    
    for i in range (64):
        for j in range(64):
            #Nos quedamos con bloques 8x8
            bloque = I_JPEG[i*8:(i*8)+8, j*8:(j*8)+8]
            
            #Hacemos la DCT a cada valor del bloque
            A = dctn(bloque, norm='ortho')
            
            #Cuantizamos los valores obtenidos por la DCT 2D
            A_quant=np.round(A/q)*q
            
            #Realizamos la transformada inversa IDCT 2D
            I_JPEG_reconstruida[i*8:(i*8)+8, j*8:(j*8)+8] = idctn(A_quant, norm='ortho') 
    
    
    
    figure, (ax1, ax2, ax3)=plt.subplots(1, 3, figsize=(15, 5))
    ax1.imshow(I_JPEG, cmap='gray')
    ax2.set_axis_off()
    ax1.set_title('Imagen original')
    ax2.imshow(I_JPEG_reconstruida, cmap='gray')
    ax2.set_axis_off()
    ax2.set_title('Imagen reconstruida, JPEG cuantizar q=' + str(q))
    #ax3.imshow(np.abs(I_JPEG-I_JPEG_reconstruida), cmap='gray')
    ax3.imshow(I_JPEG_reconstruida[176:176+128,101:101+128], cmap='gray') #Hacemos zoom en las plumas
    ax2.set_axis_off()
    ax3.set_title('Zoom para q=' + str(q))
    plt.tight_layout()
    
    PSNR_JPEG[ind]=psnr(I_JPEG,I_JPEG_reconstruida)
    ind=ind+1
    
plt.figure(ind+1)
plt.xlabel('q')
plt.ylabel('PSNR (dB)')
plt.title('PNSR en función del cuanto')
plt.plot(q0, PSNR_JPEG)
plt.grid()

#Vemos que conforme aumenta el cuanto, la calidad de la imagen sufre, ya que van quedando
#menos coeficientes sin ser 0 y nos vamos quedando con menos componentes en frecuencia

#%%CODIFICADOR JPEG: Sin cuantizador

#En primer lugar leemos la imagen y la convertimos en escala de grises
I_JPEG2 = io.imread('lena.bmp', as_gray=True).astype(np.double)*255 -128

#Si las componentes de la imagen no se pueden dividir bloques de 8, añadimos ceros
resto = len(I_JPEG2) % 8 
if resto == 0:
    pass
else:
    nceros = 8 - resto
    ceros = np.zeros(nceros)
    I_JPEG2 = np.append(I_JPEG2,ceros, axis = 0)
    I_JPEG2 = np.append(I_JPEG2, ceros, axis = 1)



#Inicializamos variables
bloque = np.zeros((8,8)) #Combinación lineal de las 64 frecuencias
I_JPEG2_reconstruida = np.zeros_like(I_JPEG2)
PSNR_JPEG2=np.zeros(10)
ind=0

#Seleccionamos el cuanto


for q in range (5, 55, 5):
    
    q0[ind]=q
    bloque = np.zeros((8,8)) #Combinación lineal de las 64 frecuencias
    
    for i in range (64):
        for j in range(64):
            #Nos quedamos con bloques 8x8
            bloque = I_JPEG2[i*8:(i*8)+8, j*8:(j*8)+8]
            
            #Hacemos la DCT a cada valor del bloque
            A2 = dctn(bloque, norm='ortho')
            
            #creamos la mascara para los valores que la cuantizacion vale 0
            
            mask=(np.round(A2/q)*q==0) #Nos crea una matriz donde si el valor
                                      #cuantizado es 0 habra un TRUE y si no
                                      #nos habra un FALSE
            
            #Creamos la matriz con los valores resultantes de hacer la DCT, a
            #excepcion de los valores que al cuantizarlos son cero que los
            #mantenemos cero
            B = A2
            B[mask] = 0
            
            ceros = busca_ceros(B)
            #print(ceros)
            
            
            #Realizamos la transformada inversa DCT 2D
            I_JPEG2_reconstruida[i*8:(i*8)+8, j*8:(j*8)+8] = idctn(B, norm='ortho') 
            
    
    figure, (ax1, ax2, ax3)=plt.subplots(1, 3, figsize=(15, 5))
    ax1.imshow(I_JPEG2, cmap='gray')
    ax2.set_axis_off()
    ax1.set_title('Imagen original')
    ax2.imshow(I_JPEG2_reconstruida, cmap='gray')
    ax2.set_axis_off()
    ax2.set_title('Imagen reconstruida, JPEG sin cuantizar q=' + str(q))
    ax3.imshow(I_JPEG2_reconstruida[176:176+128,101:101+128], cmap='gray')
    ax2.set_axis_off()
    ax3.set_title('Error')
    plt.tight_layout()
    
    PSNR_JPEG2[ind]=psnr(I_JPEG2,I_JPEG2_reconstruida)
    ind=ind+1
    
plt.figure(ind+1)
plt.xlabel('q')
plt.ylabel('PSNR (dB)')
plt.plot(q0, PSNR_JPEG2)
plt.title('PNSR en función del cuanto')
plt.grid()

#%% Representamos todos los tres casos de PSNR

mask = np.zeros_like(q0)
mask[:] = 1
PSNR=mask*PSNR
plt.figure(ind+1)
plt.xlabel('q')
plt.ylabel('PSNR (dB)')
plt.title('PSNR en funcion del cuanto (q)')
plt.plot(q0, PSNR_JPEG2)
plt.plot(q0, PSNR_JPEG)
plt.plot(q0, PSNR)
plt.legend(['Sin cuantizar', 'Cuantizando', 'Muestreo zonal'])
#Observamos cómo la PNSR del muestreo zonal se mantiene constante, porque siempre
#cogemos los mismos componentes de baja frecuencia y las PSNR de JPEG tienen una curva 
#similar, la cual se adapta a los componentes en frecuencia en cada caso
            