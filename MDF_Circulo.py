# -*- coding: utf-8 -*-
"""
Created on Fri May 15 21:42:20 2026

@author: gusta
"""

import numpy as np
import time

# CONFIGURACIÓN EXPERIMENTAL
num_ejecuciones = 10
tiempos_pre = []
tiempos_res = []

# CONFIGURACIÓN ESTATICA (Fuera para evitar re-alojamiento de RAM en cada iteración)
L = 2.0         
N = 101         
R = 0.8         
h = L / (N - 1) 
max_iter = 10000       
tolerancia = 1e-5     
epsilon_0 = 8.854e-12          
rho_valor = 1e-9               

print(f"Iniciando {num_ejecuciones} ejecuciones del MDF...\n")

for i in range(num_ejecuciones):
    
    #TIEMPO DE PRE-PROCESAMIENTO 
    start_pre = time.time()
    
    x = np.linspace(-L/2, L/2, N)
    y = np.linspace(-L/2, L/2, N)
    X, Y = np.meshgrid(x, y)
    mask = X**2 + Y**2 < R**2
    phi = np.zeros((N, N)) # Reiniciamos la solución a cero
    f = np.zeros((N, N))  
    f[mask] = rho_valor / epsilon_0 
    
    tiempos_pre.append(time.time() - start_pre)
    
    #TIEMPO DE RESOLUCIÓN
    start_res = time.time()
    
    for it in range(max_iter):
        phi_old = phi.copy()
        phi[1:-1, 1:-1] = 0.25 * (phi_old[1:-1, 0:-2] + phi_old[1:-1, 2:] + 
                                  phi_old[0:-2, 1:-1] + phi_old[2:, 1:-1] + 
                                  h**2 * f[1:-1, 1:-1])
        phi[~mask] = 0 
        error = np.max(np.abs(phi - phi_old))
        if error < tolerancia:
            break
            
    tiempos_res.append(time.time() - start_res)
    print(f"Ejecución {i+1}: Conv. en {it} it. | Res: {tiempos_res[-1]:.4f}s")

#PROMEDIOS Y MEMORIA
promedio_pre = np.mean(tiempos_pre)
promedio_res = np.mean(tiempos_res)

mem_phi_mb = phi.nbytes / (1024**2)
mem_total_mb = (phi.nbytes + f.nbytes + X.nbytes + Y.nbytes) / (1024**2)

print(f"\n--- DATOS LISTOS PARA TABLA 1 (MDF CÍRCULO) ---")
print(f"Nodos totales: {N*N}")
print(f"RAM Estática (Solo Solución Phi): {mem_phi_mb:.6f} MB")
print(f"RAM Estática (Infraestructura Total): {mem_total_mb:.6f} MB")
print(f"Tiempo Pre-procesamiento Promedio: {promedio_pre:.6f} s")
print(f"Tiempo Resolución Promedio: {promedio_res:.6f} s")
print(f"Total: {promedio_pre + promedio_res:.6f} s")