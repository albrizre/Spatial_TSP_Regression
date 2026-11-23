# A short check of the public functions, data alignment and exact small solver.
source("R/setup.R")
source("R/generation.R")
source("R/predict_tsp_optimal_length.R")
source("R/training_TSP_model.R")

square_coords <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
stopifnot(isTRUE(all.equal(exact_small_tsp(as.matrix(dist(square_coords))), 4)))

opt <- read.csv("data/tsplib_optima.csv")
coord <- read.csv("data/tsplib_coords.csv")
counts <- table(coord$instance)
stopifnot(nrow(opt) == 77, length(counts) == 77,
          !anyDuplicated(opt$instance),
          all(is.finite(opt$opt_exact) & opt$opt_exact > 0),
          all(as.integer(counts[opt$instance]) == opt$n))

# Check the TSPLIB parser and pre-fitted prediction on an actual study instance.
cc <- as.matrix(coord[coord$instance == "a280", c("x", "y")])
file <- tempfile(fileext = ".tsp")
writeLines(c("NAME : a280", "TYPE : TSP", "EDGE_WEIGHT_TYPE : EUC_2D",
             paste("DIMENSION :", nrow(cc)), "NODE_COORD_SECTION",
             paste(seq_len(nrow(cc)), cc[, 1], cc[, 2]), "EOF"), file)
stopifnot(isTRUE(all.equal(unname(read_tsp_coords(file)[, 1]), unname(cc[, 1]))))
pred <- predict_tsp_optimal_length(file)
stopifnot(nrow(pred) == 1, pred$fit > 0,
          pred$lwr <= pred$fit, pred$fit <= pred$upr)
unlink(file)

data <- read.csv("data/simulated_instances.csv")
stopifnot(nrow(data) == 150, sum(data$sd >= 1) == 148)
# Check that repeated CV produces one error per instance and repetition,
# and that the same seed reproduces the same errors.
toy <- data.frame(x = 1:24, type = rep(c("Random", "Regular"), each = 12))
toy$opt <- 10 + 2 * toy$x + sin(toy$x)
E1 <- repeated_cv(opt ~ x, toy, k = 3, reps = 4)
E2 <- repeated_cv(opt ~ x, toy, k = 3, reps = 4)
stopifnot(identical(dim(E1), c(24L, 4L)), identical(E1, E2),
          all(is.finite(E1)), all(E1 >= 0))
message("Smoke test completed.")

# Confirm the public functions match the manuscript's model and input API.
model3 <- readRDS("models/default_tsp_model.rds")
stopifnot(inherits(model3, "tsp_model"), model3$model == "scale_logCE",
          model3$n_train == 77L)
# Fit on actual archived TSPLIB optima (no heuristic lengths).
example_train <- read.csv("data/tsplib_optima.csv")
example_train$opt <- example_train$opt_exact
model_user <- make_tsp_model(example_train, model = "scale_logCE")
stopifnot(inherits(model_user, "tsp_model"), model_user$n_train == 77L,
          is.finite(model_user$smearing))
