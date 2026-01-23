# Multicomponent Diffusion in Basaltic Melts  
**Temperature-Independent Eigenvector Matrix & Multicomponent Diffusion Calculator**

## Overview

This GitHub repository contains **open-access data files and MATLAB codes** associated with the paper:

**Bai, B., Zhang, Y. (2025).**  
*Multicomponent diffusion in basaltic melts: A temperature independent eigenvector matrix, and a multicomponent diffusion calculator.*  
[**Geochimica et Cosmochimica Acta**, **407**, 133–143](https://doi.org/10.1016/j.gca.2025.09.002).


The full archival version of this dataset is hosted on the **University of Michigan Deep Blue Data Repository**.  
**Deep Blue link / DOI:** *(https://doi.org/10.7302/77zk-m861)*

---

## Scientific Background

Multicomponent diffusion in natural silicate melts is a key process in magma mixing and evolution. In an **N-component** system, multicomponent diffusion is characterized by an **(N−1) × (N−1) diffusion matrix [D]**.

- The **eigenvectors** of [D] define **(N−1) eigen-components** that diffuse independently.
- The **eigenvalues** correspond to the diffusion coefficients of these eigen-components.

Diffusion eigenvectors in **8-component basaltic melts** (SiO<sub>2</sub>–TiO<sub>2</sub>–Al<sub>2</sub>O<sub>3</sub>–FeO–MgO–CaO–Na<sub>2</sub>O–K<sub>2</sub>O) have been shown to be approximately **temperature independent**. In this study, diffusion data from multiple temperatures are simultaneously fitted using a **single eigenvector matrix [Q]**, further supporting the temperature independence of diffusion eigenvectors in basaltic melts.


---

## What’s Included

### Data Files

1. **`MultiComponentDif_calculator_v1.0.xlsx`**  
   An Excel-based calculator for computing multicomponent diffusion profiles in diffusion couples.

2. **`Bai_and_Zhang_2025_data_in_figures.xlsx`**  
   Data used to generate all figures in the associated publication.

### MATLAB Codes

3. **`code_global_fit`**  
   MATLAB program (`BFGS_main.m`) that uses **BFGS optimization** to obtain:
   - a temperature-independent eigenvector matrix [Q]
   - three sets of eigenvalues corresponding to **1260 °C, 1350 °C, and 1500 °C**  
   Includes subroutines and diffusion data from **27 diffusion couple experiments**.  
   See the `Readme` inside the folder for details.

4. **`code_fit_Z6_vs_x_in_BS13&14C`**  
   MATLAB program (`Newton_main.m`) that uses the **Newton method** to fit the **Z <sub>6</sub> vs. x** diffusion profile in experiment **BS13&14C**, assuming the diffusion coefficient **λ<sub>6</sub>** increases exponentially with **Z<sub>6</sub>**.  
   See the `Readme` inside the folder for details.


---

## Methodology Summary

1. **Eigenvector and eigenvalue fitting (global fit)**
   - MATLAB implementation using the **BFGS method**
   - Fitting data from **27 diffusion couple experiments**
   - 8-component basaltic melts with similar compositions
   - Experiments conducted at **1260 °C, 1350 °C, and 1500 °C**
   - One shared eigenvector matrix \([Q]\) is used across all temperatures

2. **Concentration-dependent diffusivity fitting (single profile)**
   - MATLAB implementation using the **Newton method** (`Newton_main.m`)
   - Diffusivity is assumed to depend exponentially on concentration

---

## Requirements

- MATLAB (recent versions recommended)
- Optimization Toolbox (recommended)

---

## Related Publication

- **Bai, B., Zhang, Y. (2025).**  
  *Multicomponent diffusion in basaltic melts: A temperature independent eigenvector matrix, and a multicomponent diffusion calculator.*  
  **Geochimica et Cosmochimica Acta**, **407**, 133–143.

---

## How to Cite

### Data citation (Deep Blue)

> Bai, B., Zhang, Y. **Data of the paper "Multicomponent diffusion in basaltic melts: A temperature independent eigenvector matrix, and a multicomponent diffusion calculator"** \[Data set\]. University of Michigan – Deep Blue Data.

**Deep Blue DOI / URL:** *(https://doi.org/10.7302/77zk-m861)*

### Code citation (this GitHub repository)

If you use or cite this repository, please include:
- Authors: Bobo Bai and Youxue Zhang  
- Repository name  
- Year  
- GitHub URL  
- Release tag or commit hash  

---

## Contact

**Bobo Bai**  
University of Michigan  
Email: bbai@umich.edu

