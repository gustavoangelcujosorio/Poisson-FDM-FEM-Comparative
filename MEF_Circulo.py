# -*- coding: utf-8 -*-
"""
Created on Thu May  7 01:09:13 2026

@author: gusta
"""

from dolfin import *
from mshr import *
import numpy as np
import time

# CONFIGURACIÓN EXPERIMENTAL
num_ejecuciones = 1
tiempos_pre = []
tiempos_res = []

print(f"Iniciando {num_ejecuciones} ejecuciones del MEF para promediado estadístico...\n")

for i in range(num_ejecuciones):
    #INICIO PRE-PROCESAMIENTO
    # Incluye generación de malla, espacios y formas variacionales
    start_pre = time.time()
    
    R = 0.8
    rho_valor = 1e-9
    epsilon_0 = 8.854e-12
    A = rho_valor / epsilon_0
    x0, y0, sigma = 0.0, 0.0, 0.05

    # 1. Malla y Espacio funcional
    circulo = Circle(Point(0, 0), R)
    malla = generate_mesh(circulo, 32) 
    V = FunctionSpace(malla, 'P', 1)

    # 2. Definición del problema variacional
    f = Expression('A * exp(-(pow(x[0]-x0, 2) + pow(x[1]-y0, 2)) / (2*pow(sigma, 2)))',
                   degree=2, A=A, x0=x0, y0=y0, sigma=sigma)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx
    
    def frontera(x, on_boundary): return on_boundary
    bc = DirichletBC(V, Constant(0.0), frontera)

    # 3. Ensamble de la matriz (Parte del pre-procesamiento)
    A_mat = assemble(a)
    b_vec = assemble(L)
    bc.apply(A_mat, b_vec)
    
    tiempos_pre.append(time.time() - start_pre)

    # INICIO RESOLUCIÓN
    # Resolución del sistema lineal ya ensamblado
    start_res = time.time()
    
    u_sol = Function(V)
    solve(A_mat, u_sol.vector(), b_vec)
    
    tiempos_res.append(time.time() - start_res)
    
    if i == 0:
        nodos_mef = malla.num_vertices()
        print(f"Nota: Malla generada con {nodos_mef} nodos (vértices).")

# CÁLCULO DE MEMORIA (IEEE 754)
# En MEF reportamos el tamaño del vector solución (Grados de Libertad)
# u_sol.vector().get_local() nos da el arreglo de floats
dofs = u_sol.vector().size()
mem_solucion_mb = (dofs * 8) / (1024**2) 

# Estimación de memoria de la Matriz de Rigidez (Formato disperso)
# Aproximadamente 7 valores no nulos por fila para elementos P1 en 2D
mem_matriz_mb = (dofs * 7 * 8) / (1024**2)

# --- CÁLCULO DE PROMEDIOS ---
promedio_pre = np.mean(tiempos_pre)
promedio_res = np.mean(tiempos_res)

print(f"\n--- RESULTADOS FINALES MEF PARA LA TABLA 1 (Promedio de {num_ejecuciones} ejecuciones) ---")
print(f"Nodos (Vértices): {malla.num_vertices()}")
print(f"Grados de Libertad (DoF): {dofs}")
print(f"Memoria Estática (Vector Solución): {mem_solucion_mb:.6f} MB")
print(f"Memoria Estática (Matriz Rigidez Est.): {mem_matriz_mb:.6f} MB")
print(f"Tiempo Pre-procesamiento Promedio: {promedio_pre:.6f} s")
print(f"Tiempo Resolución Promedio: {promedio_res:.6f} s")
print(f"Total: {promedio_pre + promedio_res:.6f} s")