# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:19:12 2026

@author: usuario
"""

def deltaE_T(i,j):
    ip=i
    im=i
    jp=j
    jm=j
    
    """
    Comprobamos la i y asignamos valores a ip (i+1) y a im (i-1)
    """
    if i==0:
        ip=i+1
        im=N-1
    elif i==N-1:
        ip=0
        im=j-1
    else:
        ip=i+1
        im=i-1
    
    """
    Comprobamos la j y asignamos valores a jp (j+1) y a jm (j-1)
    """
    if j==0:
        jp=j+1
        jm=N-1
    elif j==N-1:
        jp=0
        jm=j-1
    else:
        jp=j+1
        jm=j-1
        
    deentret=np.exp((-2*s[tempn,i,j]*(s[tempn,im,j]+s[tempn,ip,j]+s[tempn,i,jm]+s[tempn,i,jp]))/t)
    return deentret

#########################################################################################
#########################################################################################
#########################################################################################

import numpy as np

"""
Especificamos tamaño y generamos temperaturas entre 0 y 5 cada 0.5
"""
N=16
aT=np.arange(0,5,0.5)
aT[0]=0.1
print(aT)

"""
Abrimos ficheros para todo: los primeros 10 para la configuración desordenada y los siguientes 10 para la ordenada
"""
temp0=open('temp0.txt','w')
temp1=open('temp0.5.txt','w')
temp2=open('temp1.txt','w')
temp3=open('temp1.5.txt','w')
temp4=open('temp2.txt','w')
temp5=open('temp2.5.txt','w')
temp6=open('temp3.txt','w')
temp7=open('temp3.5.txt','w')
temp8=open('temp4.txt','w')
temp9=open('temp4.5.txt','w')
temp10=open('temp5.txt','w')

temp0o=open('temp0o.txt','w')
temp1o=open('temp0.5o.txt','w')
temp2o=open('temp1o.txt','w')
temp3o=open('temp1.5o.txt','w')
temp4o=open('temp2o.txt','w')
temp5o=open('temp2.5o.txt','w')
temp6o=open('temp3o.txt','w')
temp7o=open('temp3.5o.txt','w')
temp8o=open('temp4o.txt','w')
temp9o=open('temp4.5o.txt','w')
temp10o=open('temp5o.txt','w')

"""
Bucle para generar la configuración aleatoria de espines para las distintas temperaturas
"""
s=np.zeros((10,N,N))
smm=np.ones((10,N,N))           #Configuración ordenada de espín 1
prob=np.random.random((10,N,N))

for i in range(10):
    for j in range(N):
        for k in range(N):
            if prob[i,j,k] < 0.5:
                s[i,j,k]=-1
            else:
                s[i,j,k]=1
    #print(s[i,:,:])

tempn=0         #Valor que especificará en que posición de aT estamos

"""
Generamos un array con 2 valores, del cual sacaremos después el mínimo entre 1 y 
exp(-delta[E]/T)
"""
cond=np.ones(2)
condmm=np.ones(2)

"""
Creamos el bucle grande y empezamos los cálculos por cada temperatura
"""
vez=0
estacionario=np.array([False,False,False,False,False,False,False,False,False,False,False],dtype=bool)

for t in aT:
    #rand=np.random.randint(0,N,2)
    print(f"Temperatura: {t}")
    nmm=0
    desvest=0
    for i in np.arange(0,N**2*1e4):
        """
        Generamos la posición aleatoria que podrá cambiar el espin y dseta
        """
        rand=np.random.randint(0,N,2)
        dseta=np.random.random(1)
        
        """
        Comprobamos si cambia calculando el mínimo, comparamos con dseta y cambiamos si
        se cumple
        """
        cond[1]=deltaE_T(rand[0],rand[1])

        if np.min(cond)>dseta:
            s[tempn,rand[0],rand[1]]*=-1
            
        """
        Repetimos para la configuración ordenada con el objetivo de la magnetización media
        """
        rand=np.random.randint(0,N,2)
        dseta=np.random.random(1)
        condmm[1]=deltaE_T(rand[0],rand[1])
        
        if np.min(condmm)>dseta:
            smm[tempn,rand[0],rand[1]]*=-1
                
        """
        Guardamos los datos cada NxN pasos
        """
        if i%(N**2)==0:
            print(f'pmc {i/N**2}')      #Vemos cuántos pmc tenemos por seguridad
            if i>=9e4:
                nmm+=1
                
            if tempn==0:
                for j in range(N):
                    for k in range(N-1):
                        temp0.write(f'{s[tempn,j,k]}, ')
                        temp0o.write(f'{smm[tempn,j,k]}, ')
                    temp0.write(f'{s[tempn,j,N-1]}\n')
                    temp0o.write(f'{smm[tempn,j,N-1]}\n')
                temp0.write("\n")
                temp0o.write("\n")
            elif tempn==1:
                for j in range(N):
                    for k in range(N-1):
                        temp1.write(f'{s[tempn,j,k]}, ')
                        temp1o.write(f'{smm[tempn,j,k]}, ')
                    temp1.write(f'{s[tempn,j,N-1]}\n')
                    temp1o.write(f'{smm[tempn,j,N-1]}\n')
                temp1.write("\n")
                temp1o.write("\n")
            elif tempn==2:
                for j in range(N):
                    for k in range(N-1):
                        temp2.write(f'{s[tempn,j,k]}, ')
                        temp2o.write(f'{smm[tempn,j,k]}, ')
                    temp2.write(f'{s[tempn,j,N-1]}\n')
                    temp2o.write(f'{smm[tempn,j,N-1]}\n')
                temp2.write("\n")
                temp2o.write("\n")
            elif tempn==3:
                for j in range(N):
                    for k in range(N-1):
                        temp3.write(f'{s[tempn,j,k]}, ')
                        temp3o.write(f'{smm[tempn,j,k]}, ')
                    temp3.write(f'{s[tempn,j,N-1]}\n')
                    temp3o.write(f'{smm[tempn,j,N-1]}\n')
                temp3.write("\n")
                temp3o.write("\n")
            elif tempn==4:
                for j in range(N):
                    for k in range(N-1):
                        temp4.write(f'{s[tempn,j,k]}, ')
                        temp4o.write(f'{smm[tempn,j,k]}, ')
                    temp4.write(f'{s[tempn,j,N-1]}\n')
                    temp4o.write(f'{smm[tempn,j,N-1]}\n')
                temp4.write("\n")
                temp4o.write("\n")
            elif tempn==5:
                for j in range(N):
                    for k in range(N-1):
                        temp5.write(f'{s[tempn,j,k]}, ')
                        temp5o.write(f'{smm[tempn,j,k]}, ')
                    temp5.write(f'{s[tempn,j,N-1]}\n')
                    temp5o.write(f'{smm[tempn,j,N-1]}\n')
                temp5.write("\n")
                temp5o.write("\n")
            elif tempn==6:
                for j in range(N):
                    for k in range(N-1):
                        temp6.write(f'{s[tempn,j,k]}, ')
                        temp6o.write(f'{smm[tempn,j,k]}, ')
                    temp6.write(f'{s[tempn,j,N-1]}\n')
                    temp6o.write(f'{smm[tempn,j,N-1]}\n')
                temp6.write("\n")
                temp6o.write("\n")
            elif tempn==7:
                for j in range(N):
                    for k in range(N-1):
                        temp7.write(f'{s[tempn,j,k]}, ')
                        temp7o.write(f'{smm[tempn,j,k]}, ')
                    temp7.write(f'{s[tempn,j,N-1]}\n')
                    temp7o.write(f'{smm[tempn,j,N-1]}\n')
                temp7.write("\n")
                temp7o.write("\n")
            elif tempn==8:
                for j in range(N):
                    for k in range(N-1):
                        temp8.write(f'{s[tempn,j,k]}, ')
                        temp8o.write(f'{smm[tempn,j,k]}, ')
                    temp8.write(f'{s[tempn,j,N-1]}\n')
                    temp8o.write(f'{smm[tempn,j,N-1]}\n')
                temp8.write("\n")
                temp8o.write("\n")
            elif tempn==9:
                for j in range(N):
                    for k in range(N-1):
                        temp9.write(f'{s[tempn,j,k]}, ')
                        temp9o.write(f'{smm[tempn,j,k]}, ')
                    temp9.write(f'{s[tempn,j,N-1]}\n')
                    temp9o.write(f'{smm[tempn,j,N-1]}\n')
                temp9.write("\n")
                temp9o.write("\n")
            elif tempn==10:
                for j in range(N):
                    for k in range(N-1):
                        temp10.write(f'{s[tempn,j,k]}, ')
                        temp10o.write(f'{smm[tempn,j,k]}, ')
                    temp10.write(f'{s[tempn,j,N-1]}\n')
                    temp10o.write(f'{smm[tempn,j,N-1]}\n')
                temp10.write("\n")
                temp10o.write("\n")

            #print(f'pmc {pmc}')      #Vemos cuántos pmc tenemos por seguridad
    tempn+=1    #Nos movemos a la siguiente temperatura
    
"""
Cerramos los ficheros
"""
temp0.close()
temp1.close()
temp2.close()
temp3.close()
temp4.close()
temp5.close()
temp6.close()
temp7.close()
temp8.close()
temp9.close()
temp10.close()
temp0o.close()
temp1o.close()
temp2o.close()
temp3o.close()
temp4o.close()
temp5o.close()
temp6o.close()
temp7o.close()
temp8o.close()
temp9o.close()
temp10o.close()