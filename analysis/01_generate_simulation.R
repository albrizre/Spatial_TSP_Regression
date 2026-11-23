# 01. Generate the main simulated dataset (150 instances).
# Optional: the reported analyses use the archived data in data/.
# Outputs are kept apart from the archived inputs used in the paper.
source("R/setup.R")
source("R/generation.R")

folder <- "generated/simulation"
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(folder, "features.csv")))
  stop("A generated experiment already exists. Move that folder before regenerating.")

# The 5 target sizes and 10 replicates are crossed with 3 pattern types.
# The observation window is fixed to a 10 x 10 square in every case.
W <- square(10)
sizes <- c(10, 20, 50, 100, 200)
r_dist <- c(3.6, 2.4, 1.6, 1, 0.7)
grid <- expand.grid(r = 1:10, n_exp = sizes)

# Cauchy clustering, homogeneous Poisson and sequential inhibition.
# The parameters are those of the original experiment (Section 3.1.1).
# The realized sample size pp$n can differ from the target n.
simulate_pattern <- function(type, n) {
  repeat {
    pp <- switch(type,
      clust = rCauchy(kappa = 5 / 100, scale = 0.25, mu = n / 5, win = W, n.cond = n),
      hom = rpoispp(lambda = n / 100, win = W),
      reg = rSSI(r = r_dist[match(n, sizes)], win = W))
    if (pp$n >= 2) return(pp)
  }
}

features <- list()
points <- list()
for (type in c("clust", "hom", "reg")) {
  # The original scripts reset the seed separately for each spatial process.
  set.seed(12345)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    name <- paste0(type, "_ins_", i)
    pp <- simulate_pattern(type, g$n_exp)
    # Use the known observation window to compute the CE statistic.
    ce <- clarkevans.test(pp, correction = "Donnelly")
    instance <- save_pattern(pp, name, folder)
    # Preserve the order of the original random draws (heuristic before tours).
    heuristic <- solve_TSP(instance$tsp)
    s <- random_tours_sd(instance$tsp)
    features[[name]] <- data.frame(instance_name = name, count = i,
      window_size = 10, n_exp = g$n_exp, n_sim = pp$n, r = g$r,
      sd = s, opt_heuristic = tour_length(heuristic),
      R = unname(ce$statistic), p = ce$p.value, Type = type)
    points[[name]] <- data.frame(instance_name = name, instance$coords)
  }
}
# Archive both covariates and point coordinates for subsequent exact solving.
features <- do.call(rbind, features)
points <- do.call(rbind, points)
write.csv(features, file.path(folder, "features.csv"), row.names = FALSE)
write.csv(points, file.path(folder, "coordinates.csv"), row.names = FALSE)

# Obtain exact lengths. This is done after generating all point patterns
# so that it does not affect the original sequence of random draws.
concorde_ok <- concorde_available()
optima <- lapply(features$instance_name, function(name) {
  cc <- points[points$instance_name == name, c("x", "y")]
  D <- as.matrix(dist(cc))
  if (nrow(D) <= 14) {
    # Exact dynamic programming is feasible for these small instances.
    len <- exact_small_tsp(D)
    method <- "dynamic_programming"
  } else if (concorde_ok) {
    tour <- tryCatch(solve_TSP(TSP(as.dist(D)), method = "concorde"),
                     error = function(e) NULL)
    len <- if (is.null(tour)) NA_real_ else tour_length(tour)
    method <- if (is.null(tour)) "missing" else "concorde"
  } else {
    # If no exact solver is configured, leave the length missing.
    len <- NA_real_
    method <- "missing"
  }
  data.frame(instance_name = name, opt_exact = len, method = method)
})
write.csv(do.call(rbind, optima), file.path(folder, "optima.csv"), row.names = FALSE)
message("New design saved in ", folder, ". The archived data have not been changed.")
