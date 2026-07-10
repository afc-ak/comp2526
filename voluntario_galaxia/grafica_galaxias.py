import numpy as np
import matplotlib.pyplot as plt

N = 100

# Cargar el fichero
datos = np.loadtxt("estrellas_lichu.txt")

# Primera línea (t = 0)
instante0 = datos[0]

# Extraer coordenadas
x = instante0[0::2]
y = instante0[1::2]

# Representar
plt.figure(figsize=(7,7))
plt.plot(x, y, 'k.', markersize=3)

plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("")
plt.ylabel("")
plt.title("")
plt.grid(False)

plt.show()