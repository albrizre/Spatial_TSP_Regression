# Exploiting spatial structure to estimate TSP tour lengths

This repository contains the R code and data used to reproduce the main results of the paper *Exploiting spatial structure to estimate TSP tour lengths*. It also provides functions to train the proposed regression models and predict tour lengths for new Euclidean TSP instances.

## Analysis scripts

The scripts in `analysis/` are numbered and can be run individually from the repository root.

**Simulation studies**

- `01_generate_simulation.R`: Generates clustered, random and regular spatial instances.
- `02_simulation.R`: Fits and compares the regression models on the main simulated dataset.
- `03_generate_cluster_sensitivity.R`: Generates additional clustered instances for the sensitivity study.
- `04_cluster_sensitivity.R`: Evaluates the fitted models under different clustering parameters.

**TSPLIB analysis and figures**

- `05_tsplib.R`: Analyses the TSPLIB instances and compares the proposed models with alternative methods.
- `06_figures.R`: Produces the figures in the paper from the saved analysis results.

## Reproducing the results

Install the required R packages once:

```r
install.packages(c("spatstat", "TSP", "ggplot2", "ggrepel", "MASS", "openxlsx"))
```

The original datasets and optimal tour lengths are included in `data/`. To reproduce the reported analyses, run:

```r
source("analysis/02_simulation.R")
source("analysis/04_cluster_sensitivity.R")
source("analysis/05_tsplib.R")
source("analysis/06_figures.R")
```

Scripts 01 and 03 document and reproduce the simulation designs, but do not need to be rerun to analyse the original instances. Shared functions are stored in `R/`, and results and figures are saved in `results/` and `figures/`, respectively.

## Training and prediction functions

Two functions in `R/` support the application of the models to other instances:

- `training_TSP_model.R`: Fits Models 1–3 using user-provided TSP instances and their **optimal tour lengths**.
- `predict_tsp_optimal_length.R`: Predicts tour lengths and provides prediction intervals, using either a user-trained model or the supplied Model 3 fitted to the 77 TSPLIB instances.

A working example of both functions is provided in `examples/train_and_predict.R`:

```r
source("examples/train_and_predict.R")
```
