# -*- coding: utf-8 -*-
"""
Created on Sat May 9 16:36:39 2026

@author: usuario
"""

import numpy as np
import matplotlib.pyplot as plt

"""
Abrimos los ficheros necesarios
"""
pos=open('posiciones.txt','w')

"""
Definimos unas constantes y calculamos Delta y mu
"""
N=4

G=6.67e-11
MT=5.9736e24
ML=0.07349e24
m=20e6
dTL=3.844e8
omega=2.6617e-6
RT=6.37816e6
RL=1.7374e6

Delta=G*MT/(dTL**3)
mu=ML/MT

t=0

"""
Creamos un vector posición para la Luna, escogemos la h y calculamos f, vector de las funciones principales
"""
dln=np.zeros(100001)
rnave=np.zeros((int(1e5+1),2))

rluna=np.zeros((int(1e5+1),2))
r=np.zeros(2)
r[0]=1
r[1]=0
rluna[0,0]=r[0]*np.cos(r[1])
rluna[0,1]=r[0]*np.sin(r[1])

f=np.zeros(4)
fh2=np.zeros(4)
fh2aux=np.zeros(4)
ffinal=np.zeros(4)

f[0]=RT/dTL                                                  #perteneciente a r
f[1]=np.pi/4                                             #perteneciente a fi
f[2]=2*np.sqrt((2*G*MT)/RT)/dTL*np.cos(np.pi/4-f[1])         #perteneciente a pr
f[3]=f[0]*2*np.sqrt((2*G*MT)/RT)/dTL*np.sin(np.pi/4-f[1])    #perteneciente a pfi

fh2[0]=RT/dTL                                                  #perteneciente a r
fh2[1]=np.pi/4                                             #perteneciente a fi
fh2[2]=2*np.sqrt((2*G*MT)/RT)/dTL*np.cos(np.pi/4-fh2[1])         #perteneciente a pr
fh2[3]=fh2[0]*2*np.sqrt((2*G*MT)/RT)/dTL*np.sin(np.pi/4-fh2[1])    #perteneciente a pfi
fh2aux=np.copy(fh2)

ffinal[0]=RT/dTL                                                  #perteneciente a r
ffinal[1]=np.pi/4                                             #perteneciente a fi
ffinal[2]=2*np.sqrt((2*G*MT)/RT)/dTL*np.cos(np.pi/4-ffinal[1])         #perteneciente a pr
ffinal[3]=f[0]*2*np.sqrt((2*G*MT)/RT)/dTL*np.sin(np.pi/4-ffinal[1])    #perteneciente a pfi

rnave[0,0]=f[0]*np.cos(f[1])
rnave[0,1]=f[0]*np.sin(f[1])
dln[0]=np.sqrt(f[0]**2 +r[0]**2 -2*f[0]*r[0]*np.cos(f[1]-r[1]))

h=59.987654321   #h en segundos, cercana a 1 minuto
errormax=100   ####poner a 100

"""
Creamos los vectores k sup n, cada uno con cuatro coordenadas para cada función e yn, vector de las derivadas de f
"""
k1=np.zeros(N)
k2=np.zeros(N)
k3=np.zeros(N)
k4=np.zeros(N)
k1h2=np.zeros(N)
k2h2=np.zeros(N)
k3h2=np.zeros(N)
k4h2=np.zeros(N)

yn=np.zeros(N)
ynh2=np.zeros(N)
ynaux=np.zeros(N)
ynh2aux=np.zeros(N)

yn[0]=2*np.sqrt((2*G*MT)/RT)/dTL
yn[1]=0
yn[2]=((f[0]*yn[1])**2)/((f[0]/dTL)**3)-Delta*(1/(f[0]**2)+(mu/((np.sqrt(1+f[0]**2-2*f[0]*np.cos(f[1]-omega*t)))**3))*(f[0]-np.cos(f[1]-omega*t)))
yn[3]=-((Delta*mu*f[0])/((np.sqrt(1+f[0]**2-2*f[0]*np.cos(f[1]-omega*t)))**3))*np.sin(f[1]-omega*t)

ynh2[0]=2*np.sqrt((2*G*MT)/RT)/dTL
ynh2[1]=0
ynh2[2]=((f[0]*yn[1])**2)/((f[0]/dTL)**3)-Delta*(1/(f[0]**2)+(mu/((np.sqrt(1+f[0]**2-2*f[0]*np.cos(f[1]-omega*t)))**3))*(f[0]-np.cos(f[1]-omega*t)))
ynh2[3]=-((Delta*mu*f[0])/((np.sqrt(1+f[0]**2-2*f[0]*np.cos(f[1]-omega*t)))**3))*np.sin(f[1]-omega*t)

ynaux=np.copy(yn)
ynh2aux=np.copy(ynh2)

#plt.plot(f[0]*np.cos(f[1]),f[0]*np.sin(f[1]))
"""
Empezamos el bucle
"""
s=0
epsilon=0
epsilonn=np.zeros(N)

#print(t)
for j in range(int(1e5)):
    t+=h
    #print(t)
    
    """
    Establecemos las condiciones iniciales del paso con el resultado final
    """
    f=np.copy(ffinal)
    
    """
    Calculamos las k sup x sub n; x,n E [0,3] para h
    """
    #Índice 0: r·
    #Índice 1: fi·
    #Índice 2: pr·
    #Índice 3: pfi·
    
    k1[0]=h*f[2]
    k1[1]=h*f[3]/((f[0])**2)
    k1[2]=h*(f[3]**2/(f[0]**3)-Delta*(1/(f[0]**2)+mu/(f[0]**3)*(f[0]-np.cos(f[1]-omega*t))))
    k1[3]=h*(Delta*mu*f[0]/((np.sqrt(1+f[0]**2 -2*f[0]*np.cos(f[1]-omega*t)))**3)*np.sin(f[1]-omega*t))
    
    k2[0]=h*(f[2]+k1[2]/2)
    k2[1]=h*((f[3]+k1[3]/2)/((f[0]+k1[0]/2)**2))
    k2[2]=h*((f[3]+k1[3]/2)**2/((f[0]+k1[0]/2)**3)-Delta*(1/((f[0]+k1[0]/2)**2)+mu/((f[0]+k1[0]/2)**3)*((f[0]+k1[0]/2)-np.cos(f[1]+k1[1]/2-omega*t))))
    k2[3]=h*(Delta*mu*(f[0]+k1[0]/2)/((np.sqrt(1+(f[0]+k1[0]/2)**2 -2*(f[0]+k1[0]/2)*np.cos(f[1]+k1[1]/2-omega*t)))**3))*np.sin(f[1]+k1[1]/2-omega*t)
    
    k3[0]=h*(f[2]+k2[2]/2)
    k3[1]=h*((f[3]+k2[3]/2)/((f[0]+k2[0]/2)**2))
    k3[2]=h*((f[3]+k2[3]/2)**2/((f[0]+k2[0]/2)**3)-Delta*(1/((f[0]+k2[0]/2)**2)+mu/((f[0]+k2[0]/2)**3)*((f[0]+k2[0]/2)-np.cos(f[1]+k2[1]/2-omega*t))))
    k3[3]=h*(Delta*mu*(f[0]+k2[0]/2)/((np.sqrt(1+(f[0]+k2[0]/2)**2 -2*(f[0]+k2[0]/2)*np.cos(f[1]+k2[1]/2-omega*t)))**3))*np.sin(f[1]+k2[1]/2-omega*t)
    
    k4[0]=h*(f[2]+k3[2])
    k4[1]=h*((f[3]+k3[3])/((f[0]+k3[0])**2))
    k4[2]=h*((f[3]+k3[3])**2/((f[0]+k3[0])**3)-Delta*(1/((f[0]+k3[0])**2)+mu/((f[0]+k3[0])**3)*((f[0]+k3[0])-np.cos(f[1]+k3[1]-omega*t))))
    k4[3]=h*(Delta*mu*(f[0]+k3[0])/((np.sqrt(1+(f[0]+k3[0])**2 -2*(f[0]+k3[0])*np.cos(f[1]+k3[1]-omega*t)))**3))*np.sin(f[1]+k3[1]-omega*t)
    
    """
    Calculamos yn para paso h
    """
    yn[0]+=(1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
    yn[1]+=(1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    yn[2]+=(1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    yn[3]+=(1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
    
    for p in range(2):
        """
        Calculamos las k sup x sub n; x,n E [0,3] para paso h/2
        """
        #Índice 0: r·
        #Índice 1: fi·
        #Índice 2: pr·
        #Índice 3: pfi·
        
        k1h2[0]=h/2*fh2[2]
        k1h2[1]=h/2*fh2[3]/((fh2[0])**2)
        k1h2[2]=h/2*(fh2[3]**2/(fh2[0]**3)-Delta*(1/(fh2[0]**2)+mu/(fh2[0]**3)*(fh2[0]-np.cos(fh2[1]-omega*t))))
        k1h2[3]=h/2*(Delta*mu*fh2[0]/((np.sqrt(1+fh2[0]**2 -2*fh2[0]*np.cos(fh2[1]-omega*t)))**3)*np.sin(fh2[1]-omega*t))
        
        k2h2[0]=h/2*(fh2[2]+k1h2[2]/2)
        k2h2[1]=h/2*((fh2[3]+k1h2[3]/2)/((fh2[0]+k1h2[0]/2)**2))
        k2h2[2]=h/2*((fh2[3]+k1h2[3]/2)**2/((fh2[0]+k1h2[0]/2)**3)-Delta*(1/((fh2[0]+k1h2[0]/2)**2)+mu/((fh2[0]+k1h2[0]/2)**3)*((fh2[0]+k1h2[0]/2)-np.cos(fh2[1]+k1h2[1]/2-omega*t))))
        k2h2[3]=h/2*(Delta*mu*(fh2[0]+k1h2[0]/2)/((np.sqrt(1+(fh2[0]+k1h2[0]/2)**2 -2*(fh2[0]+k1h2[0]/2)*np.cos(fh2[1]+k1h2[1]/2-omega*t)))**3))*np.sin(fh2[1]+k1h2[1]/2-omega*t)
        
        k3h2[0]=h/2*(fh2[2]+k2h2[2]/2)
        k3h2[1]=h/2*((fh2[3]+k2h2[3]/2)/((fh2[0]+k2h2[0]/2)**2))
        k3h2[2]=h/2*((fh2[3]+k2h2[3]/2)**2/((fh2[0]+k2h2[0]/2)**3)-Delta*(1/((fh2[0]+k2h2[0]/2)**2)+mu/((fh2[0]+k2h2[0]/2)**3)*((fh2[0]+k2h2[0]/2)-np.cos(fh2[1]+k2h2[1]/2-omega*t))))
        k3h2[3]=h/2*(Delta*mu*(fh2[0]+k2h2[0]/2)/((np.sqrt(1+(fh2[0]+k2h2[0]/2)**2 -2*(fh2[0]+k2h2[0]/2)*np.cos(fh2[1]+k2h2[1]/2-omega*t)))**3))*np.sin(fh2[1]+k2h2[1]/2-omega*t)
        
        k4h2[0]=h/2*(fh2[2]+k3h2[2])
        k4h2[1]=h/2*((fh2[3]+k3h2[3])/((fh2[0]+k3h2[0])**2))
        k4h2[2]=h/2*((fh2[3]+k3h2[3])**2/((fh2[0]+k3h2[0])**3)-Delta*(1/((fh2[0]+k3h2[0])**2)+mu/((fh2[0]+k3h2[0])**3)*((fh2[0]+k3h2[0])-np.cos(fh2[1]+k3h2[1]-omega*t))))
        k4h2[3]=h/2*(Delta*mu*(fh2[0]+k3h2[0])/((np.sqrt(1+(fh2[0]+k3h2[0])**2 -2*(fh2[0]+k3h2[0])*np.cos(fh2[1]+k3h2[1]-omega*t)))**3))*np.sin(fh2[1]+k3h2[1]-omega*t)
        
        """
        Calculamos yn para paso h/2
        """
        ynh2[0]+=(1/6)*(k1h2[0]+2*k2h2[0]+2*k3h2[0]+k4h2[0])
        ynh2[1]+=(1/6)*(k1h2[1]+2*k2h2[1]+2*k3h2[1]+k4h2[1])
        ynh2[2]+=(1/6)*(k1h2[2]+2*k2h2[2]+2*k3h2[2]+k4h2[2])
        ynh2[3]+=(1/6)*(k1h2[3]+2*k2h2[3]+2*k3h2[3]+k4h2[3])
        
        """
        Recalculamos f para h/2
        """
        fh2[0]+=(m/MT)*ynh2[0]*t
        fh2[1]+=ynh2[1]*t
        fh2[2]+=ynh2[2]*t
        fh2[3]+=ynh2[3]*t
    
    """
    Calculamos el error para cada ecuación y nos quedamos con el mayor
    """
    epsilonn=16*np.abs(ynh2-yn)/15
    epsilon=np.max(epsilonn)
    
    """
    Seguimos con el algoritmo empezando el bucle para cambiar h
    """
    s=np.max(((epsilon/errormax)**2,1e-8))
    hmax=h/s
    
    yn=np.copy(ynaux)
    ynh2=np.copy(ynh2aux)
    
    while s>2:
        fh2=np.copy(fh2aux)
        h=hmax
        
        """
        Calculamos las k sup x sub n; x,n E [1,4] para paso h
        """
        #Índice 0: r·
        #Índice 1: fi·
        #Índice 2: pr·
        #Índice 3: pfi·
        
        k1[0]=h*f[2]
        k1[1]=h*f[3]/((f[0])**2)
        k1[2]=h*(f[3]**2/(f[0]**3)-Delta*(1/(f[0]**2)+mu/(f[0]**3)*(f[0]-np.cos(f[1]-omega*t))))
        k1[3]=h*(Delta*mu*f[0]/((np.sqrt(1+f[0]**2 -2*f[0]*np.cos(f[1]-omega*t)))**3)*np.sin(f[1]-omega*t))
        
        k2[0]=h*(f[2]+k1[2]/2)
        k2[1]=h*((f[3]+k1[3]/2)/((f[0]+k1[0]/2)**2))
        k2[2]=h*((f[3]+k1[3]/2)**2/((f[0]+k1[0]/2)**3)-Delta*(1/((f[0]+k1[0]/2)**2)+mu/((f[0]+k1[0]/2)**3)*((f[0]+k1[0]/2)-np.cos(f[1]+k1[1]/2-omega*t))))
        k2[3]=h*(Delta*mu*(f[0]+k1[0]/2)/((np.sqrt(1+(f[0]+k1[0]/2)**2 -2*(f[0]+k1[0]/2)*np.cos(f[1]+k1[1]/2-omega*t)))**3))*np.sin(f[1]+k1[1]/2-omega*t)
        
        k3[0]=h*(f[2]+k2[2]/2)
        k3[1]=h*((f[3]+k2[3]/2)/((f[0]+k2[0]/2)**2))
        k3[2]=h*((f[3]+k2[3]/2)**2/((f[0]+k2[0]/2)**3)-Delta*(1/((f[0]+k2[0]/2)**2)+mu/((f[0]+k2[0]/2)**3)*((f[0]+k2[0]/2)-np.cos(f[1]+k2[1]/2-omega*t))))
        k3[3]=h*(Delta*mu*(f[0]+k2[0]/2)/((np.sqrt(1+(f[0]+k2[0]/2)**2 -2*(f[0]+k2[0]/2)*np.cos(f[1]+k2[1]/2-omega*t)))**3))*np.sin(f[1]+k2[1]/2-omega*t)
        
        k4[0]=h*(f[2]+k3[2])
        k4[1]=h*((f[3]+k3[3])/((f[0]+k3[0])**2))
        k4[2]=h*((f[3]+k3[3])**2/((f[0]+k3[0])**3)-Delta*(1/((f[0]+k3[0])**2)+mu/((f[0]+k3[0])**3)*((f[0]+k3[0])-np.cos(f[1]+k3[1]-omega*t))))
        k4[3]=h*(Delta*mu*(f[0]+k3[0])/((np.sqrt(1+(f[0]+k3[0])**2 -2*(f[0]+k3[0])*np.cos(f[1]+k3[1]-omega*t)))**3))*np.sin(f[1]+k3[1]-omega*t)
        
        """
        Calculamos yn para paso h
        """
        yn[0]+=(1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        yn[1]+=(1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        yn[2]+=(1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
        yn[3]+=(1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
        
        for p in range(2):
            """
            Calculamos las k sup x sub n; x,n E [1,4] para paso h/2
            """
            #Índice 0: r·
            #Índice 1: fi·
            #Índice 2: pr·
            #Índice 3: pfi·
            
            k1h2[0]=h/2*fh2[2]
            k1h2[1]=h/2*fh2[3]/((fh2[0])**2)
            k1h2[2]=h/2*(fh2[3]**2/(fh2[0]**3)-Delta*(1/(fh2[0]**2)+mu/(fh2[0]**3)*(fh2[0]-np.cos(fh2[1]-omega*t))))
            k1h2[3]=h/2*(Delta*mu*fh2[0]/((np.sqrt(1+fh2[0]**2 -2*fh2[0]*np.cos(fh2[1]-omega*t)))**3)*np.sin(fh2[1]-omega*t))
            
            k2h2[0]=h/2*(fh2[2]+k1h2[2]/2)
            k2h2[1]=h/2*((fh2[3]+k1h2[3]/2)/((fh2[0]+k1h2[0]/2)**2))
            k2h2[2]=h/2*((fh2[3]+k1h2[3]/2)**2/((fh2[0]+k1h2[0]/2)**3)-Delta*(1/((fh2[0]+k1h2[0]/2)**2)+mu/((fh2[0]+k1h2[0]/2)**3)*((fh2[0]+k1h2[0]/2)-np.cos(fh2[1]+k1h2[1]/2-omega*t))))
            k2h2[3]=h/2*(Delta*mu*(fh2[0]+k1h2[0]/2)/((np.sqrt(1+(fh2[0]+k1h2[0]/2)**2 -2*(fh2[0]+k1h2[0]/2)*np.cos(fh2[1]+k1h2[1]/2-omega*t)))**3))*np.sin(fh2[1]+k1h2[1]/2-omega*t)
            
            k3h2[0]=h/2*(fh2[2]+k2h2[2]/2)
            k3h2[1]=h/2*((fh2[3]+k2h2[3]/2)/((fh2[0]+k2h2[0]/2)**2))
            k3h2[2]=h/2*((fh2[3]+k2h2[3]/2)**2/((fh2[0]+k2h2[0]/2)**3)-Delta*(1/((fh2[0]+k2h2[0]/2)**2)+mu/((fh2[0]+k2h2[0]/2)**3)*((fh2[0]+k2h2[0]/2)-np.cos(fh2[1]+k2h2[1]/2-omega*t))))
            k3h2[3]=h/2*(Delta*mu*(fh2[0]+k2h2[0]/2)/((np.sqrt(1+(fh2[0]+k2h2[0]/2)**2 -2*(fh2[0]+k2h2[0]/2)*np.cos(fh2[1]+k2h2[1]/2-omega*t)))**3))*np.sin(fh2[1]+k2h2[1]/2-omega*t)
            
            k4h2[0]=h/2*(fh2[2]+k3h2[2])
            k4h2[1]=h/2*((fh2[3]+k3h2[3])/((fh2[0]+k3h2[0])**2))
            k4h2[2]=h/2*((fh2[3]+k3h2[3])**2/((fh2[0]+k3h2[0])**3)-Delta*(1/((fh2[0]+k3h2[0])**2)+mu/((fh2[0]+k3h2[0])**3)*((fh2[0]+k3h2[0])-np.cos(fh2[1]+k3h2[1]-omega*t))))
            k4h2[3]=h/2*(Delta*mu*(fh2[0]+k3h2[0])/((np.sqrt(1+(fh2[0]+k3h2[0])**2 -2*(fh2[0]+k3h2[0])*np.cos(fh2[1]+k3h2[1]-omega*t)))**3))*np.sin(fh2[1]+k3h2[1]-omega*t)
            
            """
            Calculamos yn para paso h/2
            """
            ynh2[0]+=(1/6)*(k1h2[0]+2*k2h2[0]+2*k3h2[0]+k4h2[0])
            ynh2[1]+=(1/6)*(k1h2[1]+2*k2h2[1]+2*k3h2[1]+k4h2[1])
            ynh2[2]+=(1/6)*(k1h2[2]+2*k2h2[2]+2*k3h2[2]+k4h2[2])
            ynh2[3]+=(1/6)*(k1h2[3]+2*k2h2[3]+2*k3h2[3]+k4h2[3])
            
            """
            Recalculamos f para h/2
            """
            fh2[0]+=(m/MT)*ynh2[0]*t
            fh2[1]+=ynh2[1]*t
            fh2[2]+=ynh2[2]*t
            fh2[3]+=ynh2[3]*t
        
        """
        Calculamos el error para cada ecuación y nos quedamos con el mayor
        """
        epsilonn=16*np.abs(ynh2-yn)/15
        epsilon=np.max(epsilonn)
        
        """
        Seguimos con el algoritmo
        """
        s=np.max(((epsilon/errormax)**2,1e-8))
        hmax=h/s
        
        yn=np.copy(ynaux)
        ynh2=np.copy(ynh2aux)
    
    ffinal=np.copy(fh2)
    fh2aux=np.copy(fh2)
    ynaux=np.copy(yn)
    ynh2aux=np.copy(ynh2)
    
    if h<hmax:
        h*=2
    
    ###############################################################################
    
    rnave[j+1,0]=ffinal[0]*np.cos(ffinal[1])
    rnave[j+1,1]=ffinal[0]*np.sin(ffinal[1])
    
    r[1]=omega*t
    rluna[j+1,0]=r[0]*np.cos(r[1])
    rluna[j+1,1]=r[0]*np.sin(r[1])
    
    dln[j+1]=np.sqrt(ffinal[0]**2 +r[0]**2 -2*ffinal[0]*r[0]*np.cos(ffinal[1]-r[1]))

    pos.write(f"{rluna[j][0]}\t{rluna[j][1]}\t{rnave[j][0]}\t{rnave[j][1]}\t{dln[j]}\n")

#plt.plot(range(int(100001)),dln)
#plt.show()

pos.close()
#plt.plot(rluna[:1000,0],rluna[:1000,1],'.')    
#plt.plot(rnave[:1000,0],rnave[:1000,1])
   
#plt.plot(rnave[:500,0],rnave[:500,1],'.')
    
#plt.plot(rnave[:100,0],rnave[:100,1],'.')