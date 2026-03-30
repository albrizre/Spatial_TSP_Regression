This repository contains the R codes that allow reproducing the main results of the paper *Exploiting the spatial structure to predict optimal travelling salesperson tour lengths*.

The scripts are organized sequentially into two main groups: **Simulation Studies** and **TSPLIB Analysis**.

### 1. Simulation Studies (1a – 1e)
These R scripts allow conducting the simulation studies described in the paper.

* **1a Sim study clust.R**: Generates clustered instances with a Neyman-Scott process.
* **1b Sim study hom.R**: Generates random instances with a Poisson homogeneous process.
* **1c Sim study reg.R**: Generates regular instances with a simple sequential inhibition process.
* **1d Linear regression analysis.R**: Performs a linear regression analysis on optimal tour lengths, using the standard deviation as a covariate (with or without interaction with the Clark-Evans statistic).
* **1e Functional regression analysis.R**: Performs a functional regression analysis on optimal tour lengths, using the cumulative distribution function as a functional covariate (with or without interaction with the Clark-Evans statistic).

### 2. TSPLIB Analysis (2a – 2b)
These R scripts allow appling the proposed methods to the TSPLIB benchmark library.

* **2a TSPlib data generation.R**: Reads TSPLIB files and obtains node coordinates, converting them into a format suitable for spatial analysis in R (e.g., `ppp` objects).
* **2b TSPlib analysis.R**: Performs the classification of the TSPLIB instances based on their spatial structure and a linear regression analysis on optimal tour lengths, using the standard deviation as a covariate (with or without interaction with the Clark-Evans statistic).

Finally, we provide two extra R functions that would allow the user to implement the main methods presented with different sets of instances:

* **training_TSP_model.R**: This function allows training the linear regression models presented in the paper with a set of selected instances.
* **predict_tsp_optimal_length.R**: This function allows predicting the optimal TSP tour length for a given set of instances. The user can specify the coefficients of the linear regression model from a model fitted with function training_TSP_model.R. Otherwise, the models fitted with the TSPlib library are employed for the prediction.
