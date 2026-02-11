This repository contains the R code used associated with the paper...

The scripts are organized sequentially into two main groups: **Simulation Studies** and **TSPLIB Empirical Analysis**.

### 1. Simulation Studies (1a – 1f)
These scripts establish the baseline by simulating various spatial point processes and applying functional analysis.

* **1a Sim study Poisson clust.R**: 
* **1b Sim study Poisson hom.R**: 
* **1c Sim study Poisson inhib.R**:
* **1d Analysis all types.R**: 
* **1e Analysis all types FDA.R**: 
* **1f FPCA.R**: 


### 2. TSPLIB Analysis (2a – 2b)
These scripts apply the proposed methods to the TSPLIB benchmark library.

* **2a TSPlib data generation.R**: Parses raw TSPLIB files and extracts node coordinates, converting them into a format suitable for spatial analysis in R (e.g., `ppp` objects).
* **2b TSPlib analysis.R**: The core empirical script. It applies the classification tests (including the Clark-Evans test) to the TSPLIB instances and generates the results discussed in the paper.
