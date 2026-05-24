# A Comparative Analysis Between FDM and FEM: Toward a Mesh-Free Alternative

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![FEniCS](https://img.shields.io/badge/FEniCS-supported-orange.svg)](https://fenicsproject.org/)

This repository contains the official source code and supplementary material for the article: **"Una comparativa entre MDF y MEF: hacia una alternativa \textit{mesh-free}"**.

The project provides the complete numerical modeling of the electrostatic potential governed by the Poisson equation using two traditional mesh-based approaches: the Finite Difference Method (FDM) and the Finite Element Method (FEM). The code focuses on benchmarking computational cost (pre-processing vs. resolution time), static RAM consumption (Float64 precision), and geometric truncation errors near curved boundaries.

## 📁 Repository Structure

The experimental setup is divided into the following Python scripts:

* `MDF_Placa.py`: FDM implementation for an orthogonal 2D domain (plate) using Jacobi relaxation. Computes baseline metrics.
* `MDF_Circulo.py`: FDM implementation mapped to a curved boundary (circle). Executes a statistical average of 10 runs to evaluate the computational load of the *staircase effect*.
* `MEF_Circulo.py`: FEM implementation using unstructured Delaunay meshes via the **FEniCS** framework. Separates the cost of stiffness matrix assembly (pre-processing) from the linear system resolution.
* `Visualizacion_MDF.py`: Rendering script that plots the analytical mask and internal/boundary nodes to evidence geometric truncation.
* `Visualizacion_MEF.py`: Rendering script that plots the conforming unstructured mesh using `mshr`.

## ⚙️ Requirements and Dependencies

To reproduce the experiments, a Python 3.x environment is required. The standard scientific packages used are:
* `numpy`
* `matplotlib`
* `sys`
* `time`

**Note on FEM:** The Finite Element Method scripts strictly require the **FEniCS** computing platform and the `mshr` module. You can find installation instructions on the [official FEniCS website](https://fenicsproject.org/download/).

## 🚀 How to Run

1. Clone this repository:
```bash
   git clone [https://github.com/gustavoangelcujosorio/Poisson-FDM-FEM-Comparative.git](https://github.com/gustavoangelcujosorio/Poisson-FDM-FEM-Comparative.git)
   cd Poisson-FDM-FEM-Comparative
