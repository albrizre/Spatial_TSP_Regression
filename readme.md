This repository contains the R codes that allows reproducing the main results of the paper *Predicting optimal TSP tour lengths using spatial statistics and functional data methods*, by Álvaro Briz-Redón, Teresa León, and Juanjo Peiró.

The scripts are organized sequentially into two main groups: **Simulation Studies** and **TSPLIB Analysis**.

### 1. Simulation Studies (1a – 1f)
These R scripts establish the baseline by simulating various spatial point processes and applying functional data analysis.

* **1a Sim study clust.R**: Generates clustered instances with a Neyman-Scott process.
* **1b Sim study hom.R**: Generates random instances with a Poisson homogeneous process.
* **1c Sim study reg.R**: Generates regular instances with a simple sequential inhibition process.
* **1d Linear regression analysis.R**: Performs a linear regression analysis on optimal tour lengths, using the standard deviation as a covariate (with or without interaction with the Clark-Evans statistic).
* **1e Functional regression analysis.R**: Performs a functional regression analysis on optimal tour lengths, using the cumulative distribution function as a functional covariate (with or without interaction with the Clark-Evans statistic).

### 2. TSPLIB Analysis (2a – 2b)
These R scripts apply the proposed methods to the TSPLIB benchmark library.

* **2a TSPlib data generation.R**: Reads TSPLIB files and obtains node coordinates, converting them into a format suitable for spatial analysis in R (e.g., `ppp` objects).
* **2b TSPlib analysis.R**: Performs the classification of the TSPLIB instances based on their spatial structure and a linear regression analysis on optimal tour lengths, using the standard deviation as a covariate (with or without interaction with the Clark-Evans statistic).  
