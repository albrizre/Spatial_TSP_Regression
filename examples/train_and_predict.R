# Example of the two public functions, using small TSPLIB instances included
# in data/tsplib_coords.csv. Run from the repository root.
source("R/training_TSP_model.R")
source("R/predict_tsp_optimal_length.R")

coords <- read.csv("data/tsplib_coords.csv")
optima <- read.csv("data/tsplib_optima.csv")
train_names <- c("berlin52", "eil51", "eil76", "eil101",
                 "ch130", "ch150", "kroA100", "kroB100")
test_name <- "a280"

# Write a few standard .tsp files to a temporary directory so that the
# example can be run without downloading anything else.
tmp_dir <- tempfile("tsp_example_")
dir.create(tmp_dir)
names_to_write <- c(train_names, test_name)
files <- file.path(tmp_dir, paste0(names_to_write, ".tsp"))
for (i in seq_along(names_to_write)) {
  name <- names_to_write[i]
  xy <- coords[coords$instance == name, c("x", "y")]
  writeLines(c(paste("NAME :", name),
               "TYPE : TSP", "EDGE_WEIGHT_TYPE : EUC_2D",
               paste("DIMENSION :", nrow(xy)), "NODE_COORD_SECTION",
               paste(seq_len(nrow(xy)), xy$x, xy$y), "EOF"), files[i])
}

# 1. Train Model 3 using known optimal lengths (never heuristic estimates).
opt_train <- optima$opt_exact[match(train_names, optima$instance)]
my_model <- training_TSP_model(files[1:length(train_names)], opt = opt_train)
print(my_model)

# 2. Predict a new instance with the fitted model, and with the bundled
# Model 3 trained on all 77 study instances. Both provide 95% intervals.
pred_user <- predict_tsp_optimal_length(files[length(files)], model = my_model)
pred_default <- predict_tsp_optimal_length(files[length(files)])
print(pred_user[, c("fit_mean", "lwr", "upr", "extrapolation")])
print(pred_default[, c("fit_mean", "lwr", "upr", "extrapolation")])

# fit_mean: smearing-corrected point prediction (used in the manuscript's
# error calculations). fit: exponentiated fitted log length, without correction.
# lwr/upr: prediction interval on the original length scale.
unlink(tmp_dir, recursive = TRUE)
