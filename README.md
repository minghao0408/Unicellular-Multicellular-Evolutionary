
# Unicellular-Multicellular Evolutionary Branching Reproduction

## Project Overview

This repository contains a comprehensive reproduction and analysis of the computational model described in the paper:

> **"Unicellular–multicellular evolutionary branching driven by resource limitations"** > *Adriano Bonforti & Ricard Solé (2022), Journal of The Royal Society Interface*

### Scientific Context
The transition from unicellular to multicellular life is a major evolutionary milestone. This project implements a **Cellular Potts Model (CPM)** coupled with evolutionary dynamics to investigate how **spatial resource limitations** alone can drive the evolution of cell-cell adhesion.

**Key Findings Reproduced:**
* **Critical Transition:** A sharp phase transition occurs at a critical resource patch size $\phi_c \approx 60$.
* **Multicellularity (MC):** For small resource patches ($\phi < \phi_c$), cells evolve negative adhesion values ($J < 0$), forming tight aggregates to stay near nutrients.
* **Unicellularity (UC):** For large resource patches ($\phi > \phi_c$), cells evolve positive adhesion values ($J > 0$), favoring dispersal.
* **Evolutionary Branching:** Complex branching dynamics and transient phenotype coexistence occur near the critical point.

---

## Repository Structure

* **`main.c`**: The core simulation engine written in C. It handles the physics (Metropolis algorithm), metabolism, and evolutionary logic. Refactored for clarity and efficiency.
* **`paper_figure_reproduction.m`**: MATLAB script to analyze data and generate the mean adhesion evolution plots (reproducing Figure 2 & 3 of the original paper).
* **`paper_figure_reproduction2.m`**: MATLAB script to visualize evolutionary branching trees and spatial snapshots (reproducing Figure 4 & 5).

---

## Getting Started

### Prerequisites
* **C Compiler**: GCC (GNU Compiler Collection) or any standard C compiler.
* **MATLAB**: For data visualization (optional; raw data is text-based).

### 1. Compilation
Navigate to the repository folder in your terminal and compile the C code using GCC:

```bash
gcc main.c -o main -lm -Wall
````

  * `-lm`: Links the math library (required for `exp`, `pow`, etc.).
  * `-Wall`: Enables all compiler warnings.

### 2\. Running the Simulation

Execute the program with the following command-line arguments:

```bash
./main [rows] [cols] [duration] [savestep] [start_J] [pop_pct] [mutation] [thresh]
```

**Standard Parameter Example (Critical Regime $\phi \approx 74$):**

```bash
./main 200 200 200000 1000 0 100 1 5
```

**Arguments Explained:**

1.  `200`: Lattice rows ($L$).
2.  `200`: Lattice columns ($L$).
3.  `200000`: Total simulation time steps ($T_{max}$).
4.  `1000`: Data saving interval (steps).
5.  `0`: Initial adhesion value ($J_0$).
6.  `100`: Initial population percentage inside the patch.
7.  `1`: Mutation enabled (1=True, 0=False).
8.  `5`: Reproduction threshold ($B_{max}$ scaled, see Implementation Note below).

> **Note on Resource Size:** To change the resource patch diameter ($\phi$), modify the `startingvarfoodpatch` variable directly in `main.c` (around line 138) and recompile.

### 3\. Visualization

After running the simulation, two data files will be generated in the root directory:

  * `outputfile.txt`: Contains spatial matrices (0=empty, 1=cell).
  * `adhspecfile.txt`: Contains adhesion values ($J$) for each cell position.

Run the provided MATLAB scripts (`.m` files) to process these files and visualize the results.

-----

## ⚠️ Important Implementation Note

### Parameter Rescaling: Code vs. Paper

Users should be aware of a numerical scaling difference between the theoretical framework described in the original paper and the actual C code implementation.

To ensure numerical stability and computational efficiency, the code employs an **Energy Rescaling Factor of 100**. This means internal energy values in the simulation are 100x larger than the theoretical units in the paper.

| Parameter | Paper Value (Theory) | Code Value (Implementation) | Notes |
| :--- | :--- | :--- | :--- |
| **Resource Consumption ($\xi_R$)** | 0.001 | **0.01** | `xi_par` in code |
| **Min. Biomass ($B_{min}$)** | 0.01 | **1.0** | Scaled by 100 |
| **Reproduction Threshold ($B_{max}$)** | 0.05 | **5.0** | Scaled by 100 |
| **Mutation Std Dev ($\sigma_m$)** | 0.005 | **0.1** | Increased to accelerate evolution |

**Conclusion:** This rescaling is linear and preserves the qualitative evolutionary dynamics described in the report.

This project was completed as part of the *Physical Methods for Biology* course (2025).

```
```
