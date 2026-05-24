# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 01:09:45 2026

@author: gusta
"""
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt

# Parametros de escala consistentes
R = 0.8

# Generar malla en FEniCS segun Listing 8
domain = Circle(Point(0, 0), R)
# Resolucion 32 genera una densidad similar a la reportada en Tabla 1
mesh = generate_mesh(domain, 32) 

print(f"Nodos reales en malla MEF: {mesh.num_vertices()}")

plt.figure(figsize=(7, 7))
plot(mesh, title=f"Malla MEF (R={R})\nAdaptación conforme mediante Triangulación de Delaunay")
plt.axis('equal')
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.show()