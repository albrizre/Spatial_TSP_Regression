# Functions used when generating new simulated instances.
source("R/features.R")
library(TSP)

# Same random-tour calculation as the original simulation.
# Each permutation is a closed cycle over all n nodes.
random_tours_sd <- function(tsp, n_random = 1000) {
  D <- as.matrix(tsp)
  n <- nrow(D)
  sd(replicate(n_random, {
    p <- sample.int(n)
    sum(D[cbind(p[-n], p[-1])]) + D[p[n], p[1]]
  }))
}

# Held-Karp dynamic programming, used for small Euclidean instances.
# Fix the first node; each state records the best path to a final node
# visiting exactly the subset represented by its bit mask.
exact_small_tsp <- function(D) {
  n <- nrow(D)
  stopifnot(n >= 2, n <= 14, ncol(D) == n)
  m <- n - 1L
  bits <- bitwShiftL(1L, seq_len(m) - 1L)
  full <- sum(bits)
  dp <- matrix(Inf, full + 1L, m)
  dp[cbind(bits + 1L, seq_len(m))] <- D[1, -1]  # single-node paths
  for (mask in seq_len(full)) {
    members <- which(bitwAnd(mask, bits) != 0L)
    if (length(members) <= 1L) next
    for (j in members) {
      previous <- mask - bits[j]
      k <- setdiff(members, j)
      dp[mask + 1L, j] <- min(dp[previous + 1L, k] + D[k + 1L, j + 1L])
    }
  }
  min(dp[full + 1L, ] + D[-1, 1])  # return from the last node to the start
}

# Check whether the configured Concorde executable can solve a small problem.
concorde_available <- function() {
  tryCatch({
    solve_TSP(TSP(dist(matrix(runif(20), 10))), method = "concorde")
    TRUE
  }, error = function(e) FALSE)
}

# Export distances with six decimal places for external exact solvers.
save_pattern <- function(pp, name, folder) {
  coords <- cbind(x = pp$x, y = pp$y)
  tsp <- TSP(dist(coords))
  write_TSPLIB(tsp, file = file.path(folder, paste0(name, ".tsp")), precision = 6)
  list(tsp = tsp, coords = coords)
}
