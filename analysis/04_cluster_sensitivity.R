# 04. Out-of-design validation on clustered instances (Section 3.1.3).
# Train on the main simulation and predict a distinct cluster design.
source("R/setup.R")
source("R/simulation_data.R")
orig <- load_simulation()

# The archived results from the revision are preferred to newly simulated data.
if (file.exists("data/cluster_sensitivity.csv")) {
  sens <- read.csv("data/cluster_sensitivity.csv")
} else if (file.exists("data/cluster_sensitivity.rda")) {
  e <- new.env()
  load("data/cluster_sensitivity.rda", envir = e)
  sens <- e$sens
} else {
  f <- "generated/cluster_sensitivity/features.csv"
  if (!file.exists(f))
    stop("Supply the archived sensitivity data or run analysis/03_generate_cluster_sensitivity.R.")
  sens <- read.csv(f)

  # Optional external EXACT lengths for these same generated patterns.
  # The CSV must contain instance_name and opt_exact (original length units).
  ext_file <- "generated/cluster_sensitivity/optima.csv"
  if (file.exists(ext_file)) {
    ext <- read.csv(ext_file)
    if (anyDuplicated(ext$instance_name)) stop("Duplicated external instance names.")
    opt <- ext$opt_exact[match(sens$instance_name, ext$instance_name)]
    ok <- is.finite(opt) & opt > 0
    sens$opt[ok] <- opt[ok]
    sens$method[ok] <- "concorde"
  }
}

# Check the design size and that all reference lengths are exact.
stopifnot(nrow(sens) == 360, !anyDuplicated(sens$instance_name))
if (any(sens$method != "concorde") || any(!is.finite(sens$opt) | sens$opt <= 0))
  stop("The sensitivity study requires exact lengths for all 360 instances.")

# Estimate the three models using only the original 148 simulated instances.
models <- list(SD = log(opt) ~ log(sd),
               SD_CE = log(opt) ~ log(sd) * R,
               n_logCE = log(opt) ~ log(n_sim) + log(R))
fits <- lapply(models, lm, data = orig)
smear <- sapply(fits, function(f) mean(exp(residuals(f))))
# Descriptive gap between construction heuristics and the exact solution.
gap <- 100 * (orig$opt_heuristic / orig$opt - 1)
print(round(t(sapply(split(gap, orig$Type), summary)), 2))

# Predict the new instances without refitting or estimating a new smearing factor.
sens <- sens[sens$sd >= 1, ]
for (m in names(fits))
  sens[[paste0("ape_", m)]] <- ape(exp(predict(fits[[m]], sens)) * smear[m], sens$opt)
sens$detected_clustered <- sens$p < 0.05 & sens$R < 1
# Summarise CE, detection rate and prediction error over the 12 cluster scenarios.
# Each scenario pools the three sizes and their replicates.
tab <- aggregate(cbind(R, detected_clustered, ape_SD, ape_SD_CE, ape_n_logCE) ~
                   parents + scale, data = sens, FUN = mean)
print(round(tab, 3))
print(round(aggregate(cbind(ape_SD, ape_SD_CE, ape_n_logCE) ~ n, sens, mean), 3))
print(round(colMeans(sens[, c("ape_SD", "ape_SD_CE", "ape_n_logCE")]), 3))
write.csv(tab, "results/cluster_sensitivity_mape.csv", row.names = FALSE)

# Use paired APEs to compare models within each cluster scenario.
# Positive differences mean that the second model is more accurate.
paired <- function(a, b) {
  groups <- split(sens, list(sens$parents, sens$scale))
  tab <- t(sapply(groups, function(d) {
    x <- d[[a]] - d[[b]]
    c(mean_diff = mean(x), share_better = mean(x > 0),
      p_wilcoxon = wilcox.test(x, exact = FALSE)$p.value)
  }))
  data.frame(model1 = a, model2 = b, scenario = rownames(tab), tab, row.names = NULL)
}
paired_tab <- rbind(paired("ape_SD", "ape_SD_CE"),
                    paired("ape_SD_CE", "ape_n_logCE"))
print(paired_tab)
write.csv(paired_tab, "results/cluster_sensitivity_paired.csv", row.names = FALSE)
# Separately refit on the sensitivity design to check coefficient stability.
# These fits are not used in the preceding prediction assessment.
refitted <- lapply(models, lm, data = sens)
print(coefficient_table(refitted))
write.csv(coefficient_table(refitted), "results/cluster_sensitivity_coefficients.csv",
          row.names = FALSE)
write.csv(sens, "results/cluster_sensitivity_predictions.csv", row.names = FALSE)
