# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:54:01 2026

@author: usuario
"""

import numpy as np

"""
Abrimos ficheros
"""
lam0=open('lam0.txt','w')
lam1=open('lam1.txt','w')
lam2=open('lam2.txt','w')
modtotal0=open('modtotal0.txt','w')
modtotal1=open('modtotal1.txt','w')
modtotal2=open('modtotal2.txt','w')

"""
Fijamos los valores iniciales y creamos los vectores temporal y espacial a recorrer, además del vector 
donde guardaremos los módulos
"""
N=500
n=N/4 #<=N/4
lambd=np.zeros(3)
for i in range(3):
    lambd[i]=0.25*(i+1)
#print(lambd)

t=np.arange(1,5e3)
modulos=np.zeros(N)

"""
Hallamos s tilda, k0 tilda, V_j tilda, fi_j0 y alfa, casi todas dentro de un bucle para
ir de lambda en lambda junto al vector A de posiciones -, 0 y +
"""
k0=(2*np.pi*n/N)
s=0.25*(k0**(-2))

alfa=np.zeros(N,dtype=complex)
gamma=np.zeros(N,dtype=complex) 
fi=np.zeros(N,dtype=complex)
chi=np.zeros(N+1,dtype=complex)
V=np.zeros(N,dtype=complex)
A=np.zeros(N,dtype=complex)
beta=np.zeros(N,dtype=complex)

for l in range(3):  #Bucle principal para pasar por todas las lambdas
    print(f'Lambda: {lambd[l]}')
    """
    Inicializamos todo lo que necesitaremos
    """
    alfa[:]=0+0j
    beta[:]=0+0j
    for i in range(int(N*2/5),int(N*3/5 +1)):
        V[i]=lambd[l]*k0**2
    A[:]=-2+((0+2*1j)/s)-V[:]
    
    for c in range(N-1,0,-1):
        gamma[c]=1/(A[c]+alfa[c])
        alfa[c-1]=-gamma[c]
    gamma[0]=1/(A[0]+alfa[0])
    
    """
    Guardamos en el límite de la caja
    """
    if l==0:
        lam0.write(f'{0}, {0}, {0}\n')
    if l==1:
        lam1.write(f'{0}, {0}, {0}\n')
    if l==2:
        lam2.write(f'{0}, {0}, {0}\n')
    normatotal=0
    
    for q in range(N):
        fi[q]=np.exp((0+1j)*(q*k0))*np.exp(((-8+0*1j)*((4*q+0*1j)-(N+0*1j))**2)/((N+0*1j)**2))
        modulos[q]=np.real(fi[q])**2+np.imag(fi[q])**2
        """
        Guardamos los resultados
        """
        if l==0:
            lam0.write(f'{np.real(fi[q])}, {np.imag(fi[q])}, {modulos[q]}\n')
        if l==1:
            lam1.write(f'{np.real(fi[q])}, {np.imag(fi[q])}, {modulos[q]}\n')
        if l==2:
            lam2.write(f'{np.real(fi[q])}, {np.imag(fi[q])}, {modulos[q]}\n')
    normatotal=np.sum(modulos)
    if l==0:
        lam0.write(f'{0}, {0}, {0}\n\n')
        modtotal0.write(f'{0}\t{normatotal}\n')
    if l==1:
        lam1.write(f'{0}, {0}, {0}\n\n')
        modtotal1.write(f'{0}\t{normatotal}\n')
    if l==2:
        lam2.write(f'{0}, {0}, {0}\n\n')
        modtotal2.write(f'{0}\t{normatotal}\n')
    #plt.plot(xreal,modulos)
    
    """
    Empezamos con el bucle de tiempos
    """
    for p in t:
        """
        Calculamos beta
        """
        for c in range(N-1,0,-1):
            beta[c-1]=gamma[c]*((((0+4*1j)*fi[c])/s)-beta[c])
        
        """
        Guardamos en el límite de la caja
        """
        if l==0:
            lam0.write(f'{0}, {0}, {0}\n')
        if l==1:
            lam1.write(f'{0}, {0}, {0}\n')
        if l==2:
            lam2.write(f'{0}, {0}, {0}\n')

        """
        Calculamos chi, posteriormente calculamos la nueva fi y el módulo
        """
        chi[:] = 0+0j
        for q in range(N):
            chi[q+1]=chi[q]*alfa[q]+beta[q]
        for q in range(N):
            fi[q]=chi[q]-fi[q]
        norm = np.sqrt(np.sum(np.abs(fi)**2))
        if norm != 0 and np.isfinite(norm):
            fi /= norm
        for q in range(N):
            modulos[q]=np.real(fi[q])**2+np.imag(fi[q])**2
            #plt.plot(xreal,modulos)
            """
            Volvemos a guardar los resultados
            """
            if l==0:
                lam0.write(f'{np.real(fi[q])}, {np.imag(fi[q])}, {modulos[q]}\n')
            if l==1:
                lam1.write(f'{np.real(fi[q])}, {np.imag(fi[q])}, {modulos[q]}\n')
            if l==2:
                lam2.write(f'{np.real(fi[q])}, {np.imag(fi[q])}, {modulos[q]}\n')
        normatotal=np.sum(modulos)
        if l==0:
            lam0.write(f'{0}, {0}, {0}\n\n')
            modtotal0.write(f'{p}\t{normatotal}\n')
        if l==1:
            lam1.write(f'{0}, {0}, {0}\n\n')
            modtotal1.write(f'{p}\t{normatotal}\n')
        if l==2:
            lam2.write(f'{0}, {0}, {0}\n\n')
            modtotal2.write(f'{p}\t{normatotal}\n')
"""
Cerramos los ficheros
"""
lam0.close()
lam1.close()
lam2.close()
modtotal0.close()
modtotal1.close()
modtotal2.close()