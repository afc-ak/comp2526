# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:53:59 2024

@author: usuario
"""

def update_plot0(frame):
    plt.cla()  # Clear current plot
    plt.imshow(data0, cmap='viridis', interpolation='nearest', vmin=0, vmax=10)
def update_plot1(frame):
    plt.cla()  # Clear current plot
    plt.imshow(data1, cmap='viridis', interpolation='nearest', vmin=0, vmax=10)
def update_plot2(frame):
    plt.cla()  # Clear current plot
    plt.imshow(data2, cmap='viridis', interpolation='nearest', vmin=0, vmax=10)
        
############################################################################################
############################################################################################
############################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N=16
#usar np.reshape para cada línea del fichero
datos0=np.loadtxt('temp0.txt')
datos1=np.loadtxt('temp2.5.txt')
datos2=np.loadtxt('temp5.txt')
data0 = np.zeros((datos0.shape[0],N,N))
data1 = np.zeros((datos1.shape[0],N,N))
data2 = np.zeros((datos2.shape[0],N,N))
data0=np.reshape(datos0, (datos0.shape[0],N,N))
data1=np.reshape(datos1, (datos1.shape[0],N,N))
data2=np.reshape(datos2, (datos2.shape[0],N,N))

#print(data0[0])

# Create a new figure and axis
fig0, ax = plt.subplots()
fig1, ax = plt.subplots()
fig2, ax = plt.subplots()
 
# Create the animation
ani0 = animation.FuncAnimation(fig0, update_plot0, frames=1000, interval=20)
ani1 = animation.FuncAnimation(fig1, update_plot1, frames=1000, interval=20)
ani2 = animation.FuncAnimation(fig2, update_plot2, frames=1000, interval=20)

# Save the animation as a GIF file
ani0.save('t0.gif', writer='pillow')
ani1.save('t2.5.gif', writer='pillow')
ani2.save('t5.gif', writer='pillow')

# Show the animation (optional)
#plt.show()