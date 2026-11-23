# 03. Generate the 360-instance sensitivity design (Section 3.1.3).
# Optional: the original 360 instances are already archived in data/.
source("R/setup.R")
source("R/generation.R")
folder <- "generated/cluster_sensitivity"
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(folder, "features.csv")))
  stop("A generated experiment already exists. Move that folder before regenerating.")

# 3 parent intensities, 4 cluster scales, 3 sizes and 10 replicates.
W <- square(10)
grid <- expand.grid(parents = c(2, 5, 10), scale = c(0.1, 0.25, 0.5, 1),
                    n = c(50, 100, 200), r = 1:10)

# Simulate extra offspring and thin to exactly n points.
# Increase the multiplier if a realization contains fewer than n points.
sim_cluster <- function(parents, scale, n) {
  mult <- 2
  repeat {
    pp <- suppressMessages(suppressWarnings(
      rCauchy(kappa = parents / 100, scale = scale,
              mu = mult * n / parents, win = W)))
    if (pp$n >= n)
      return(list(pp = pp[sample.int(pp$n, n)], how = paste0("thinned_x", mult)))
    mult <- mult * 2
  }
}

concorde_ok <- concorde_available()
set.seed(2026)
features <- list()
points <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  sim <- sim_cluster(g$parents, g$scale, g$n)
  pp <- sim$pp
  ce <- clarkevans.test(pp, correction = "Donnelly")
  name <- sprintf("clust_sens_%04d", i)
  instance <- save_pattern(pp, name, folder)
  s <- random_tours_sd(instance$tsp)

  # Retain the order of the revision script, including the heuristic
  # fallback because the solver can affect subsequent random draws.
  tour <- if (concorde_ok)
    tryCatch(solve_TSP(instance$tsp, method = "concorde"),
             error = function(e) NULL) else NULL
  if (is.null(tour)) {
    # Store a provisional heuristic result only; script 04 will reject it
    # unless an exact solution is supplied for this same realization.
    tour <- solve_TSP(instance$tsp, method = "arbitrary_insertion",
                      two_opt = TRUE, control = list(rep = 20))
    method <- "heuristic"
  } else method <- "concorde"

  features[[i]] <- data.frame(instance_name = name, g, n_sim = pp$n,
    sd = s, opt = tour_length(tour), method = method, simulation = sim$how,
    R = unname(ce$statistic), p = ce$p.value)
  points[[i]] <- data.frame(instance_name = name, instance$coords)
  if (i %% 10 == 0) message(i, "/", nrow(grid))
}
write.csv(do.call(rbind, features), file.path(folder, "features.csv"), row.names = FALSE)
write.csv(do.call(rbind, points), file.path(folder, "coordinates.csv"), row.names = FALSE)
message("New design saved in ", folder, ". Solve exported .tsp files if exact lengths are missing.")
