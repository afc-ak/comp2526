def update_plot(frame):
    plt.cla()  # Clear current plot
    for i in range(N):
        plt.xlim(-7, 7)
        plt.ylim(-7, 7)
        plt.plot(data[frame, i, 0], data[frame, i, 1], 'o', color='black', markersize=1)  # Plot current position
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axis('off')  # Turn off axes
        
############################################################################################
############################################################################################
############################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N=100

datos=np.loadtxt('estrellas_lichu.txt')
data = np.zeros((datos.shape[0],N,2))
data[:,:,0]=datos[:,0::2]
data[:,:,1]=datos[:,1::2]

# Create a new figure and axis
fig, ax = plt.subplots()
 
# Create the animation
ani = animation.FuncAnimation(fig, update_plot, frames=int(datos.shape[0]*0.4), interval=250)

# Save the animation as a GIF file
ani.save('lichu_movement.gif', writer='pillow')

# Show the animation (optional)
plt.show()
#for i in range(10):
    #plt.plot(data[:,i,0],data[:,i,1],'.')