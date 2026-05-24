# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:14:26 2026

@author: gusta
"""

import numpy as np
import time
import sys
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#INICIO PRE-PROCESAMIENTO

start_pre = time.time()

#  1. Definición de Parámetros Físicos y de Malla 
epsilon_0 = 8.854e-12
Lx, Ly = 10.0, 10.0
h = 0.1  # Tamaño del paso (Delta)

Nx = int(Lx/h) + 1
Ny = int(Ly/h) + 1
x_num = np.linspace(0, Lx, Nx)
y_num = np.linspace(0, Ly, Ny)

# 2. Métricas de Gasto Computacional 
nodos_totales = Nx * Ny
# Memoria de la matriz principal en MegaBytes
memoria_mb = (sys.getsizeof(np.zeros((Ny, Nx)) * 8)) / (1024**2) 

# 3. Inicialización y Condiciones de Frontera
Phi_num = np.zeros((Ny, Nx))
Phi_num[0, :] = 0    # Borde inferior
Phi_num[-1, :] = 0   # Borde superior
Phi_num[:, 0] = 0    # Borde izquierdo
Phi_num[:, -1] = 0   # Borde derecho

def rho(x, y):
    return 1e-9  # en nanocoulombs para valores físicos realistas

end_pre = time.time()
tiempo_pre = end_pre - start_pre

#FIN PRE-PROCESAMIENTO


# Impresión de datos estructurales de la malla (Fuera del cronómetro)
print(f"Resolución de malla: {Nx}x{Ny}")
print(f"Nodos totales (Grados de libertad): {nodos_totales}")
print(f"Memoria estimada para la matriz de potencial: {memoria_mb:.4f} MB")



#INICIO RESOLUCIÓN

start_res = time.time()

# 4. Solución Iterativa (Método de Jacobi con Criterio de Parada) 
tol = 1e-5       # Tolerancia de error
error = 1.0      # Error inicial
it = 0           # Contador de iteraciones
max_it = 10000   # Límite de seguridad

while error > tol and it < max_it:
    Phi_old = Phi_num.copy()
    
    # Bucle de actualización (Relajación de Jacobi)
    for i in range(1, Ny-1):
        for j in range(1, Nx-1):
            # Ecuación de Poisson discretizada por Diferencias Centrales
            Phi_num[i, j] = 0.25 * (Phi_old[i+1, j] + Phi_old[i-1, j] + 
                                     Phi_old[i, j+1] + Phi_old[i, j-1] + 
                                     (h**2 * rho(x_num[j], y_num[i]) / epsilon_0))
    
    # Calcular el error (Norma infinita de la diferencia)
    error = np.max(np.abs(Phi_num - Phi_old))
    it += 1

end_res = time.time()
tiempo_res = end_res - start_res

#FIN RESOLUCIÓN 


# Impresión de métricas finales
print(f"Convergencia alcanzada en {it} iteraciones.")
print(f"\n--- METRICAS DE TIEMPO REALES ---")
print(f"Tiempo de Pre-procesamiento: {tiempo_pre:.6f} segundos.")
print(f"Tiempo de Resolución:        {tiempo_res:.6f} segundos.")
print(f"Tiempo Total:                {tiempo_pre + tiempo_res:.6f} segundos.")


# Visualización de Resultados
X, Y = np.meshgrid(x_num, y_num)
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Phi_num, cmap='viridis')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('Potencial (V)')
plt.colorbar(surf)
plt.show()