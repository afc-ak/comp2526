# -*- coding: utf-8 -*-
"""
Created on Sat Apr  18 17:10:14 2026

@author: usuario
"""

def update_plot(frame):
    plt.cla()  # Clear current plot
    for i in range(N):
        plt.plot(data[:frame+1, i, 0], 
                 data[:frame+1, i, 1], 
                 color=plt.cm.viridis(i/N), alpha=0.3, linewidth=2)  # Plot previous positions with shading
        plt.plot(data[frame, i, 0], data[frame, i, 1], 'o', color=plt.cm.viridis(i/N), markersize=8)  # Plot current position
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axis('off')  # Turn off axes
        
############################################################################################
############################################################################################
############################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N=10

datos=np.loadtxt('planetas_posiciones.txt')
data = np.zeros((datos.shape[0],N,2))
data[:,:,0]=datos[:,0::2]
data[:,:,1]=datos[:,1::2]

# Create a new figure and axis
fig, ax = plt.subplots()
 
# Create the animation
ani = animation.FuncAnimation(fig, update_plot, frames=int(datos.shape[0]*0.4), interval=1)

# Save the animation as a GIF file
ani.save('planet_movement.gif', writer='pillow')

# Show the animation (optional)
plt.show()
#for i in range(10):
    #plt.plot(data[:,i,0],data[:,i,1],'.')
