## Overview

This program uses the **Newton method** to fit **Z₆ vs. x diffusion profiles**, assuming that the diffusion coefficient **λ** increases exponentially with **Z**, according to:

λ = λ₀ · exp(a · Z)

---

## How to Use the Program

This folder contains two MATLAB files:

1. **`Newton_main.m`** — main program  
2. **`Solver3.m`** — subroutine  

The main program, **`Newton_main.m`**, reads the **Z₆ vs. x** diffusion profile from the Excel file **`BS13&14C_Z6.xlsx`**.

The program assumes that the **right far-field concentration** (x → +∞) is **higher than the left far-field concentration** (x → −∞).  
If the user applies this program to a diffusion profile where the **left far-field concentration is higher**, the spatial coordinate **x must be multiplied by −1** before fitting.

---

## Fitting Parameters

The model contains **five fitting parameters**, stored in the array **`Par`**:

1. **`Par(1)`** — pre-exponential factor **λ₀** in λ = λ₀ · exp(a · Z)
2. **`Par(2)`** — exponential coefficient **a**
3. **`Par(3)`** — Matano interface position
4. **`Par(4)`** — left far-field concentration
5. **`Par(5)`** — right far-field concentration

---

## Output

After the main program finishes:

- The optimized fitting parameters are stored in the variable **`Par`**.
- The fitted **Z₆ vs. x** profile is written to the Excel file  
  **`fit_Z6_dependent.xlsx`**.

---

## Notes

For additional implementation details, numerical methods, and assumptions, please refer to the **comments and notes within each MATLAB code file**.
