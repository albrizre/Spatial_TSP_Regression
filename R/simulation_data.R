# The archived main simulation is kept separate from newly generated patterns.
load_simulation <- function(include_excluded = FALSE) {
  data <- read.csv("data/simulated_instances.csv")
  optima <- read.csv("data/simulation_optima.csv")
  # Match by instance ID rather than relying on the order of the solver output.
  if (anyDuplicated(optima$instance_name)) stop("Duplicated instance names in the optima file.")
  data$opt <- optima$opt_exact[match(data$instance_name, optima$instance_name)]
  data$opt_source <- ifelse(is.finite(data$opt) & data$opt > 0, "exact", "missing")
  if (!include_excluded) {
    # Exclude the two clustered cases with nearly constant random-tour lengths.
    # All remaining reference lengths must be exact, not construction heuristics.
    data <- data[data$sd >= 1, ]
    if ("epsilon" %in% names(data)) data <- data[data$epsilon == 0, ]
    if (any(data$opt_source != "exact"))
      stop("Missing exact reference lengths for: ",
           paste(data$instance_name[data$opt_source != "exact"], collapse = ", "),
           ". Complete data/simulation_optima.csv before analysing the simulations.")
    if (nrow(data) != 148) stop("Expected 148 retained simulated instances.")
  }
  # The simulated class is known from the generating process, not the CE test.
  data$type <- factor(c(clust = "Clustered", hom = "Random", reg = "Regular",
                        inhib = "Regular")[data$Type])
  data$log_opt <- log(data$opt)
  data$log_sd <- log(data$sd)
  data$log_n <- log(data$n_sim)
  rownames(data) <- NULL
  data
}

# Optional utility: only needed if replacing the archived exact lengths.
# Import the final exact reference lengths without modifying the simulated data.
# Excel/CSV: first column instance_name, second column the length.
# Set scale = 1e6 for the original Concorde export, 1 for unscaled lengths.
import_simulation_optima <- function(file, scale = 1e6) {
  stopifnot(is.finite(scale), scale > 0)
  ext <- if (grepl("\\.xlsx$", file, ignore.case = TRUE))
    openxlsx::read.xlsx(file) else read.csv(file)
  ext <- ext[, 1:2]
  names(ext) <- c("instance_name", "opt_exact")
  ext$opt_exact <- as.numeric(ext$opt_exact) / scale
  if (anyDuplicated(ext$instance_name)) stop("Duplicated instance names.")
  data <- read.csv("data/simulated_instances.csv")
  current <- read.csv("data/simulation_optima.csv")
  i <- match(current$instance_name, ext$instance_name)
  ok <- is.finite(ext$opt_exact[i]) & ext$opt_exact[i] > 0
  current$opt_exact[ok] <- ext$opt_exact[i[ok]]
  retained <- data$instance_name[data$sd >= 1]
  lengths <- current$opt_exact[match(retained, current$instance_name)]
  if (any(!is.finite(lengths) | lengths <= 0))
    stop("The imported file does not supply all missing reference lengths.")
  write.csv(current, "data/simulation_optima.csv", row.names = FALSE)
  message("Exact lengths saved for all 148 retained instances.")
}
