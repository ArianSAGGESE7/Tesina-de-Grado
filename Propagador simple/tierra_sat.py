import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd

# Función para leer posiciones desde un archivo CSV
def read_positions(filename):
    return pd.read_csv(filename, header=None).values

# Función para leer posiciones desde un archivo CSV y seleccionar columnas específicas
def read_positions2(filename):
    data = pd.read_csv(filename, header=None)
    return data.iloc[:, 1:4].values  # Selecciona las columnas 2 a 4 (índices 1 a 3 en Python)

# Leer las posiciones de los satélites desde los archivos CSV
satellite1_positions = read_positions('propGPS.csv')
satellite2_positions = read_positions2('estimated_positions.csv')

# Crear la figura y el eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Dibujar la Tierra
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 6371 * np.outer(np.cos(u), np.sin(v))
y = 6371 * np.outer(np.sin(u), np.sin(v))
z = 6371 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='b', alpha=0.3)

# Dibujar las posiciones de los satélites
ax.scatter(satellite1_positions[:, 0], satellite1_positions[:, 1], satellite1_positions[:, 2], color='r', label='Satélite 1')
ax.scatter(satellite2_positions[:, 0], satellite2_positions[:, 1], satellite2_positions[:, 2], color='g', label='Satélite 2')

# Etiquetas y leyenda
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

plt.show()
