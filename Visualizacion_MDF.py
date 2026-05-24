import numpy as np
import matplotlib.pyplot as plt

# Parametros exactos del Listing 5
L = 2.0
N = 101
R = 0.8
h = L / (N - 1)

x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
X, Y = np.meshgrid(x, y)

# Mascara analitica segun Listing 6
mask = X**2 + Y**2 < R**2

plt.figure(figsize=(7, 7))
# Nodos internos (los que resuelven la fisica)
plt.scatter(X[mask], Y[mask], s=1.5, c='blue', marker='o', label='Nodos internos MDF')
# Nodos de frontera (donde ocurre el escalonamiento)
plt.scatter(X[~mask], Y[~mask], s=0.5, c='gray', alpha=0.3, label='Exterior/Frontera')

# Circulo teorico para evidenciar el "staircase effect"
theta = np.linspace(0, 2*np.pi, 500)
plt.plot(R*np.cos(theta), R*np.sin(theta), 'r-', linewidth=2, label='Frontera Circular Real')

plt.title(f"Discretización MDF (L={L}, R={R}, N={N})\nEvidencia del error de escalonamiento geométrico")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.legend(loc='upper right', fontsize='small')
plt.axis('equal')
plt.grid(True, linestyle=':', alpha=0.5)
plt.show()