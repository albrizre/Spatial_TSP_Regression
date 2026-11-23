# Spatial descriptors of a Euclidean TSP instance.
library(spatstat)

# Read exactly DIMENSION rows from NODE_COORD_SECTION.
# Trim whitespace and discard the node ID (it is not a coordinate).
read_tsp_coords <- function(file) {
  lines <- trimws(readLines(file, warn = FALSE))
  header <- function(key) trimws(sub("^[^:]*:", "",
                                     grep(paste0("^", key, "\\s*:"), lines, value = TRUE)[1]))
  n <- as.integer(header("DIMENSION"))
  ewt <- toupper(header("EDGE_WEIGHT_TYPE"))
  start <- match("NODE_COORD_SECTION", lines)
  if (is.na(start) || is.na(n))
    stop("Missing NODE_COORD_SECTION or DIMENSION: ", file)
  body <- lines[seq.int(start + 1, length(lines))]
  body <- body[nzchar(body)]
  if (length(body) < n) stop("Fewer than DIMENSION nodes: ", file)
  rows <- strsplit(body[seq_len(n)], "\\s+")
  if (any(lengths(rows) != 3)) stop("Expected node ID, x and y: ", file)
  coords <- cbind(x = as.numeric(sapply(rows, "[", 2)),
                  y = as.numeric(sapply(rows, "[", 3)))
  if (any(!is.finite(coords))) stop("Invalid coordinates: ", file)
  if (!is.na(ewt) && !ewt %in% c("EUC_2D", "CEIL_2D"))
    warning("Non-Euclidean EDGE_WEIGHT_TYPE (", ewt, "): ", file)
  attr(coords, "edge_weight_type") <- ewt
  coords
}

# For TSPLIB and user-provided instances, use their bounding rectangle as W.
# Simulated patterns retain their original fixed window instead.
as_ppp <- function(coords) {
  coords <- as.matrix(coords)
  if (ncol(coords) != 2 || nrow(coords) < 2 ||
      any(!is.finite(coords)) || any(apply(coords, 2, function(x) diff(range(x))) <= 0))
    stop("A non-degenerate two-dimensional point pattern is required.")
  ppp(coords[, 1], coords[, 2], window = owin(range(coords[, 1]), range(coords[, 2])))
}

# A significant departure from CSR is classified by the direction of CE.
classify_ce <- function(R, p) {
  ifelse(p >= 0.05, "Random", ifelse(R > 1, "Regular", "Clustered"))
}

# Generate a uniformly random permutation and close the tour.
# Computing consecutive distances avoids storing the full distance matrix.
random_tour_length <- function(coords) {
  p <- sample.int(nrow(coords))
  q <- c(p, p[1])
  sum(sqrt(diff(coords[q, 1])^2 + diff(coords[q, 2])^2))
}

random_tour_sd <- function(coords, n_random = 1000) {
  sd(replicate(n_random, random_tour_length(coords)))
}

# Compute CE with Donnelly's edge correction and the area of W.
# Model 3 does not require SD, so random tours can be skipped.
tsp_features <- function(coords, with_sd = FALSE, n_random = 1000) {
  if (with_sd && (!is.numeric(n_random) || length(n_random) != 1 ||
                  !is.finite(n_random) || n_random < 2 || n_random != floor(n_random)))
    stop("n_random must be an integer of at least 2 for SD-based models.")
  pp <- as_ppp(coords)
  ce <- clarkevans.test(pp, correction = "Donnelly")
  R <- unname(ce$statistic)
  p <- ce$p.value
  out <- data.frame(n = pp$n, A = area.owin(pp$window), CE = R,
                    p_CE = p, type = tolower(classify_ce(R, p)))
  if (with_sd) out$sd <- random_tour_sd(coords, n_random)
  out
}

# The three regression models implemented by the public functions.
tsp_model_formulas <- list(
  scale_logCE = log(opt) ~ log(sqrt(n * A)) + log(CE),
  sd_ce = log(opt) ~ log(sd) * CE,
  sd = log(opt) ~ log(sd)
)

tsp_model_needs_sd <- function(model) model != "scale_logCE"

instance_features <- function(instances, model, n_random) {
  if (!length(instances)) stop("Provide at least one TSPLIB file.")
  do.call(rbind, lapply(instances, function(file) {
    cbind(instance = file, tsp_features(read_tsp_coords(file),
                 with_sd = tsp_model_needs_sd(model), n_random = n_random))
  }))
}
