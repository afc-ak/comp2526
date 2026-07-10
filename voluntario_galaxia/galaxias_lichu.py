def accel(m, r, N):
    """
    Función para calcular la aceleración tras cada paso h dados un vector de masas y la matriz 2xN de posiciones
    """
    a[:,:]=0
    for i in range(1,N):
        for j in range(N):
            if i!=j:
                a[0,i]+=-(m[j]*(r[0,i]-r[0,j]))/(np.sqrt((r[0,i]-r[0,j])**2+(r[1,i]-r[1,j])**2)**3)
                a[1,i]+=-(m[j]*(r[1,i]-r[1,j]))/(np.sqrt((r[0,i]-r[0,j])**2+(r[1,i]-r[1,j])**2)**3)
    return a

def waux(v, a, h, N):
    w=np.zeros((2,N))
    """
    Función para calcular la matriz auxiliar 2xN w tras cada paso h dados la matriz 2xN de velocidades, 
    la matriz 2xN de aceleraciones y el propio paso temporal h
    """
    for i in range(N):
        w[:,i]=v[:,i]+h/2*a[:,i]
    return w

def vel(a, w, h, N):
    v=np.zeros((2,N))
    """
    Función para calcular la aceleración tras cada paso h dados la matriz 2x10 de aceleraciones, 
    la matriz 2xN auxiliar w y el propio paso temporal h
    """
    v=np.zeros((2,N))
    for i in range(N):
        v[:,i]=w[:,i]+h/2*a[:,i]
    return v

def pos(r, v, a, h,N):
    """
    Función para calcular las posiciones tras cada paso h dados la matriz 2xN de posiciones, 
    la matriz 2x10 de velocidades, la matriz 2xN de aceleraciones y el propio paso temporal h
    """
    for i in range(N):
        r[:,i]+=h*v[:,i]+h**2/2*a[:,i]
    return r

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################

import numpy as np
import random

G=6.67392e-11               #Constante gravitación universal en N*m**2/kg**2
dSBH=2.15979e20             #Distancia Sol-Sagitario A* en m
MS=1.99e30                  #Masa solar en kg
MBH=MS*4e6                  #Masa Sagitario A* en kg
alpc=0.306601               #Conversión: 1 año luz = 0.306601 pc
mpc=3.2408e-17              #Conversión: 1 m = 3.2408*10⁽-17) pc
N=int(1e2)                  #Número de estrellas
nbrazos=int(4)              #Número de brazos que se van a representar
DSA=150e3*alpc/(dSBH*mpc)   #Diámetro de Sagitaro A*, primero en pc y después reescalado
PI=4*np.arctan(1)
a=3e3/(dSBH*mpc)
b=np.exp(np.tan(24*PI/360))
sigma_r=0.02*DSA            #Rango de desviación para las posiciones
sigma_theta=0.05            #Rango de desviaciones para los ángulos
random.seed("lichu")

posi=open('estrellas_lichu.txt','w')
periodosolw=open('periodo.txt', 'w')

"""
Inicializamos el vector masas reescalado, utilizando unos límites aproximados de las estrellas de la galaxia y 
distribuyéndolas de forma más o menos uniforme
"""
m = np.zeros(N)
m[0]=1
for i in range(1,N):
    m[i] = random.uniform(0.08, 150.0)*MS/MBH
m[N-1]=MS/MBH
#print(f"{m[1:10]}")    Para comprobar que da valores correctos

"""
Matriz posiciones reescalada: con filas acordes a coordenadas x e y de cada columna/planeta y colocamos estas posiciones en 
una matriz 2xN
"""
thetamax=np.log((DSA)/a)/b          #Ángulo máximo de desviación de la trayectoria
phi=0                               #Ángulo que determinará en qué punto de la elipse pre-rotación está la estrella
A=np.zeros(N)                       #Semieje mayor de la órbita de la estrella
e=np.zeros(N)                       #Excentricidad de las órbitas
theta=np.zeros(N)
phi=np.zeros(N)

r=np.zeros((2,N))
for k in range(1,N):
    phi[k]=random.uniform(0, 2*PI)
    if k != N-1:
        A[k]=random.uniform(a, DSA)
        e[k]=0.2 * (A[k] - a)/(DSA - a) + random.uniform(0, 0.001)      #Calculamos la excentricidad siguiendo una relación lineal
        theta[k]=np.log(A[k] / a) / b + random.uniform(-sigma_theta, sigma_theta)
        r[0,k]=(A[k] * (1 - e[k] ** 2)/(1 + e[k] * np.cos(phi[k]))) * np.cos(phi[k] + theta[k]) 
        r[1,k]=(A[k] * (1 - e[k] ** 2)/(1 + e[k] * np.cos(phi[k]))) * np.sin(phi[k] + theta[k])

    else:
        A[k]=1
        e[k]=0.085
        theta[k]=np.log(A[k] / a) / b + random.uniform(0, sigma_theta)
        r[0,k]=(A[k] * (1 - e[k] ** 2)/(1 + e[k] * np.cos(phi[k]))) * np.cos(phi[k] + theta[k]) 
        r[1,k]=(A[k] * (1 - e[k] ** 2)/(1 + e[k] * np.cos(phi[k]))) * np.sin(phi[k] + theta[k])
    r[0,k]+=(r[0,k]/DSA) * random.uniform(-sigma_r, sigma_r)
    r[1,k]+=(r[1,k]/DSA) * random.uniform(-sigma_r, sigma_r)

for i in range(N):
  posi.write(f'{r[0,i]}\t{r[1,i]}\t')
posi.write("\n")

"""
Matriz velocidad reescalada: con filas acordes a componentes v_x y v_y de cada columna/planeta

Vector auxiliar w
"""
v=np.zeros((2,N)) 
v[:,0]=0
for k in range(1,N):
    v[0,k]=-(np.sin(phi[k]) * np.cos(theta[k]) + (e[k] + np.cos(phi[k])) * np.sin(theta[k])) / (np.sqrt(A[k] * (1 - e[k] **2)))
    v[1,k]=(-np.sin(phi[k]) * np.sin(theta[k]) + (e[k] + np.cos(phi[k])) * np.cos(theta[k])) / (np.sqrt(A[k] * (1 - e[k] **2)))

w=np.zeros((2,N))

"""
Matriz aceleración: con filas acordes a componentes a_x y a_y de cada columna/planeta
"""
a=np.zeros((2,N)) 

"""
Indicamos el paso temporal h
"""
h=0.01
pasos=0
vez=0

a=accel(m, r, N)

"""
Empezamos los cálculos
"""
angulosol=np.arctan2(r[0,N-1], r[1, N-1])
angulosol_nuevo=0
delta_angulosol=0
periodosol=0
periodosol_teor=0

for t in np.arange(0, 3e3, h):
    print(f"{t}")
    """
    Recalculamos las posiciones y la matriz auxiliar w
    """
    r=pos(r, v, a, h, N)
    w=waux(v, a, h, N)
    
    angulosol_nuevo=np.arctan2(r[0,N-1], r[1, N-1])
    if angulosol_nuevo < 0:
        angulosol_nuevo += 2 * PI
    delta_angulosol=angulosol-angulosol_nuevo
    if (delta_angulosol >= -0.1) and (delta_angulosol <= 0.1) and (vez == 0):
        vez += 1
        periodosol = t * h * 1.69e10    #Calculamos el periodo en años
        periodosol_teor = 2 * PI * 1.69e10 * A[N-1] ** 1.5
        periodosolw.write(f"{periodosol}\t{periodosol_teor}")
    """
    Recalculamos la aceleración según las nuevas posiciones (posible mejora en términos de memoria)
    """
    a=accel(m, r, N)
    
    """
    Con las aceleraciones nuevas, recalculamos las velocidades
    """
    v=vel(a, w, h, N)

    """
    Aumentamos la variable para contar y guardamos las posiciones cada x pasos
    """
    pasos+=1
    if pasos%500==0:
        for i in range(N):
            posi.write(f'{r[0,i]}\t{r[1,i]}\t')
        posi.write("\n")
posi.close()
periodosolw.close()