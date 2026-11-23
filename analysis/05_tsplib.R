# 05. TSPLIB benchmark study (Section 3.2).
# Run from the root of the repository.
source("R/setup.R")
source("R/training_TSP_model.R")

# 1. TSPLIB: optimal lengths and instance descriptors
# Read and match coordinates and reference lengths by instance name.
key_of <- function(x) sub("\\.tsp$", "", x)
opt <- read.csv("data/tsplib_optima.csv")
dat <- data.frame(instance = key_of(opt$instance), opt = opt$opt_exact)
# Exclude five explicit-distance matrices and linhp318 (same nodes as lin318).
excluded <- c("bayg29", "bays29", "dantzig42", "gr120", "pa561", "linhp318")
dat <- dat[!is.na(dat$opt) & !dat$instance %in% excluded, ]
dat <- dat[order(dat$instance), ]
if (anyDuplicated(dat$instance) || any(!is.finite(dat$opt) | dat$opt <= 0))
  stop("TSPLIB instance names must be unique and optimal lengths positive.")
dat$key <- dat$instance
if (nrow(dat) != 77) stop("Expected 77 TSPLIB instances; found ", nrow(dat))
cc <- read.csv("data/tsplib_coords.csv")
# Store one coordinate matrix per instance, aligned with the length data.
coords <- lapply(split(cc[, c("x", "y")], key_of(cc$instance)), as.matrix)
missing <- setdiff(dat$instance, names(coords))
if (length(missing)) stop("Missing coordinates for: ", paste(missing, collapse = ", "))
coords <- coords[dat$instance]
reference <- opt[match(dat$instance, key_of(opt$instance)), ]

# Benchmark by Nicola et al. (2019). Use blocks to avoid an n x n matrix
# for the largest TSPLIB instances. VarP uses ordered pairs, as in the revision.
nicola_features <- function(cc, block = 1000) {
  x <- cc[, 1]
  y <- cc[, 2]
  n <- length(x)
  minP <- Inf
  maxP <- sumMin <- sumMax <- s1 <- s2 <- 0
  # Each block contains distances from a subset of nodes to all n nodes.
  # Accumulate the extrema, sums and sums of squares across blocks.
  for (start in seq(1, n, by = block)) {
    i <- start:min(n, start + block - 1)
    D <- sqrt(outer(x[i], x, "-")^2 + outer(y[i], y, "-")^2)
    D[cbind(seq_along(i), i)] <- NA  # exclude self-distances
    rmin <- apply(D, 1, min, na.rm = TRUE)
    rmax <- apply(D, 1, max, na.rm = TRUE)
    minP <- min(minP, rmin)
    maxP <- max(maxP, rmax)
    sumMin <- sumMin + sum(rmin)
    sumMax <- sumMax + sum(rmax)
    s1 <- s1 + sum(D, na.rm = TRUE)
    s2 <- s2 + sum(D^2, na.rm = TRUE)
  }
  # m counts ordered pairs; dm contains distances to the centroid.
  m <- n * (n - 1)
  dm <- sqrt((x - mean(x))^2 + (y - mean(y))^2)
  c(MinP = minP, MaxP = maxP, VarP = (s2 - s1^2 / m) / (m - 1),
    SumMinP = sumMin, SumMaxP = sumMax, MinM = min(dm), MaxM = max(dm),
    SumM = sum(dm), VarM = var(dm), VarXY = var(x) * var(y))
}

# Bounding rectangle, coordinate dispersion and midline distances,
# followed by the distance-based descriptors of Nicola et al.
geometry_features <- function(cc) {
  x <- cc[, 1]; y <- cc[, 2]
  dx <- abs(x - mean(range(x))); dy <- abs(y - mean(range(y)))
  data.frame(A = diff(range(x)) * diff(range(y)), sdx = sd(x), sdy = sd(y),
             cx = mean(abs(x - mean(x))), cy = mean(abs(y - mean(y))),
             cbar_x = mean(dx), cbar_y = mean(dy),
             cstdev_x = sd(dx), cstdev_y = sd(dy), t(nicola_features(cc)))
}

# Use archived random-tour SD by default; optionally compute 1,000 random tours.
recompute_sd <- FALSE
n_random <- 1000
set.seed(ANALYSIS_SEED)

features <- vector("list", length(coords))
for (i in seq_along(coords)) {
  cc <- coords[[i]]
  f <- tsp_features(cc, with_sd = recompute_sd, n_random = n_random)
  s <- if (recompute_sd) f$sd else reference$sd[i]
  features[[i]] <- cbind(n = f$n, sd = s, R = f$CE, p = f$p_CE,
                         geometry_features(cc))
}
dat <- cbind(dat, do.call(rbind, features))

# Verify that the raw coordinates reproduce the archived spatial descriptors.
stopifnot(all(dat$n == reference$n),
          isTRUE(all.equal(dat$A, reference$A, tolerance = 1e-10)),
          isTRUE(all.equal(dat$R, reference$CE, tolerance = 1e-6)))
if (recompute_sd && !isTRUE(all.equal(dat$sd, reference$sd, tolerance = 1e-6)))
  warning("Recomputed SD differs from the archived revision values.")
dat$type <- factor(classify_ce(dat$R, dat$p))
dat$log_opt <- log(dat$opt)
dat$log_sd <- log(dat$sd)
dat$log_n <- log(dat$n)
dat$log_A <- log(dat$A)
dat$log_sqnA <- log(sqrt(dat$n * dat$A))
dat$log_sqnS <- log(sqrt(dat$n * dat$sdx * dat$sdy))
dat$log_sqnC <- log(sqrt(dat$n * dat$cx * dat$cy))

# Terms of the Cavdar-Sokol benchmark, plus the published small-n correction.
dat$cs1 <- sqrt(dat$n * dat$cstdev_x * dat$cstdev_y)
dat$cs2 <- sqrt(dat$n * dat$sdx * dat$sdy * dat$A / (dat$cbar_x * dat$cbar_y))
correction <- ifelse(dat$n < 1000, 0.9325 * exp(0.00005298 * dat$n) -
                       0.2972 * exp(-0.01452 * dat$n), 1)
dat$cs_pub <- (2.791 * dat$cs1 + 0.2669 * dat$cs2) / correction
rownames(dat) <- NULL
print(table(dat$type))
if (!identical(as.integer(table(dat$type)), c(35L, 19L, 23L)))
  warning("TSPLIB class counts differ from the manuscript (35, 19, 23).")
write.csv(dat[, c("key", "n", "sd", "R", "p", "type", "opt")],
          "results/tsplib_covariates.csv", row.names = FALSE)

# 2. Regression models and paired repeated cross-validation
# Fit the benchmark models on the training sample of each fold.
nicola_vars <- c("n", "MinP", "MaxP", "VarP", "SumMinP", "SumMaxP",
                 "MinM", "MaxM", "SumM", "VarM", "VarXY")

# Fit on the length scale without intercept. Weights 1/L^2 target
# relative squared error and prevent very large instances dominating the fit.
fit_vars <- function(vars, d, weighted = TRUE) {
  f <- reformulate(vars, "opt", intercept = FALSE)
  if (weighted) lm(f, data = d, weights = 1 / d$opt^2) else lm(f, data = d)
}

p_values <- function(fit) {
  cf <- summary(fit)$coefficients
  setNames(cf[, 4], rownames(cf))
}

# Backward elimination: discard aliased terms and remove the least
# significant predictor until all surviving terms meet the threshold.
nicola_backward <- function(d, vars = nicola_vars, alpha = 0.05, weighted = TRUE) {
  repeat {
    fit <- fit_vars(vars, d, weighted)
    vars <- intersect(vars, names(which(!is.na(coef(fit)))))
    p <- p_values(fit)[vars]
    if (!length(vars)) stop("No estimable Nicola benchmark predictor.")
    p[is.na(p)] <- 1
    if (length(vars) == 1 || max(p) <= alpha) return(fit_vars(vars, d, weighted))
    vars <- setdiff(vars, names(which.max(p)))
  }
}

# Forward selection with a backward check after admitting each predictor.
# Both selection procedures are repeated inside every training fold.
nicola_forward <- function(d, vars = nicola_vars, alpha = 0.05, weighted = TRUE) {
  sel <- character(0)
  for (step in 1:50) {
    candidates <- setdiff(vars, sel)
    if (!length(candidates)) break
    p <- sapply(candidates, function(v) {
      pv <- p_values(fit_vars(c(sel, v), d, weighted))
      if (v %in% names(pv)) pv[[v]] else 1
    })
    p[is.na(p)] <- 1
    if (min(p) > alpha) break
    sel <- c(sel, names(which.min(p)))
    pv <- p_values(fit_vars(sel, d, weighted))[sel]
    pv[is.na(pv)] <- 1
    if (length(sel) > 1 && max(pv) > alpha)
      sel <- setdiff(sel, names(which.max(pv)))
  }
  if (!length(sel)) sel <- "SumMinP"
  fit_vars(sel, d, weighted)
}

# Cavdar and Sokol (2015), Eq. (3), re-estimated using WLS (or OLS).
cs_fit <- function(d, weighted = TRUE) fit_vars(c("cs1", "cs2"), d, weighted)

# BHH is the classical scale-only benchmark; BHH_logCE is Model 3.
# SD and SD_CE are Models 1 and 2. The remaining formulas assess
# additional covariates, interactions and class-specific effects.
models <- list(
  BHH = log_opt ~ log_sqnA,
  n_A = log_opt ~ log_n + log_A,
  CavdarSokol = cs_fit,
  SD = log_opt ~ log_sd,
  SD_n_A = log_opt ~ log_sd + log_n + log_A,
  SD_CE = log_opt ~ log_sd * R,
  SD_n_A_CE = log_opt ~ log_sd * R + log_n + log_A,
  SD_n_A_CEadd = log_opt ~ log_sd + log_n + log_A + R,
  BHH_CE = log_opt ~ log_sqnA * R,
  BHH_logCE = log_opt ~ log_sqnA + log(R),
  n_A_CE = log_opt ~ log_n + log_A + R,
  NN_sum = log_opt ~ log(SumMinP),
  NN_sum_scale = log_opt ~ log_sqnA + log(SumMinP),
  SD_class = log_opt ~ log_sd * type,
  SD_class_CE = log_opt ~ log_sd * type + log_sd * R,
  SD_classXCE = log_opt ~ log_sd * R * type
)

# Re-estimate models in each training fold. The published Cavdar-Sokol
# predictor uses fixed coefficients, so its errors are the same in every fold.
E <- list()
for (m in names(models)) {
  E[[m]] <- repeated_cv(models[[m]], dat)
}
E$CavdarSokol_pub <- matrix(ape(dat$cs_pub, dat$opt), nrow(dat), CV_REPS)
E$Nicola_back <- repeated_cv(nicola_backward, dat)
E$Nicola_fwd <- repeated_cv(nicola_forward, dat)
# Unweighted versions for comparison.
E$CavdarSokol_ols <- repeated_cv(function(d) cs_fit(d, weighted = FALSE), dat)
E$Nicola_back_ols <- repeated_cv(function(d) nicola_backward(d, weighted = FALSE), dat)
E$Nicola_fwd_ols <- repeated_cv(function(d) nicola_forward(d, weighted = FALSE), dat)

print(summary(nicola_backward(dat))$coefficients)
print(summary(nicola_forward(dat))$coefficients)
print(summary(cs_fit(dat))$coefficients)
print(summary(cs_fit(dat, weighted = FALSE))$coefficients)

# Tables of CV errors and paired comparisons (first minus second model).
pairs <- list(c("SD", "SD_CE"), c("BHH", "SD_CE"), c("n_A", "SD_CE"),
              c("CavdarSokol", "SD_CE"), c("CavdarSokol", "BHH_logCE"),
              c("CavdarSokol_pub", "BHH_logCE"), c("Nicola_back", "BHH_logCE"),
              c("Nicola_fwd", "BHH_logCE"), c("NN_sum", "BHH_logCE"),
              c("NN_sum_scale", "BHH_logCE"), c("SD_n_A", "SD_n_A_CE"),
              c("SD_CE", "SD_class_CE"), c("SD_class", "SD_class_CE"),
              c("SD_class", "SD_classXCE"), c("BHH", "BHH_CE"),
              c("BHH", "BHH_logCE"), c("n_A", "n_A_CE"),
              c("BHH_CE", "SD_n_A_CE"), c("n_A_CE", "SD_n_A_CE"),
              c("SD_n_A_CEadd", "SD_n_A_CE"), c("SD_CE", "BHH_logCE"))

tab_tsplib <- cv_table(E, dat$type)
paired_tsplib <- comparison_table(E, pairs, dat$type)
print(tab_tsplib)
print(paired_tsplib)
write.csv(tab_tsplib, "results/tsplib_cv.csv", row.names = FALSE)
write.csv(paired_tsplib, "results/tsplib_paired.csv", row.names = FALSE)

# 3. Coefficients and nested F-tests
# The F-tests compare nested log-scale regression models fitted to all instances.
models_log <- Filter(function(x) inherits(x, "formula"), models)
fits <- lapply(models_log, lm, data = dat)
print(anova(fits$SD_n_A, fits$SD_n_A_CE))
print(anova(fits$SD_CE, fits$SD_class_CE))
print(anova(fits$SD_class, fits$SD_class_CE))
print(anova(fits$BHH, fits$BHH_CE))
print(anova(fits$n_A_CE, fits$SD_n_A_CE))
print(sapply(fits, AIC))
print(sapply(fits, fit_stats))
write.csv(coefficient_table(fits), "results/tsplib_coefficients.csv", row.names = FALSE)
for (m in names(fits)) {
  cat("\n", m, "\n")
  print(round(coef_table(fits[[m]]), 4))
}
# CE times sqrt(n*A)/2 approximates the sum of nearest-neighbour distances.
print(cor(log(dat$SumMinP), log(dat$R * sqrt(dat$n * dat$A) / 2)))

# Save the fitted Model 3 separately from the pre-trained model distributed in models/.
default_tsp_model <- make_tsp_model(
  data.frame(key = dat$key, opt = dat$opt, n = dat$n, A = dat$A,
             CE = dat$R, sd = dat$sd), model = "scale_logCE")
saveRDS(default_tsp_model, "results/fitted_model3.rds")

# 4. Back-transformation bias and influence diagnostics
# Compare naive exponentiation, log-normal adjustment and Duan smearing
# under the same repeated-CV partitions.
bt <- sapply(c("naive", "lognormal", "smearing"), function(b) {
  sapply(models[c("SD", "SD_CE", "BHH_logCE")], function(m) mean(repeated_cv(m, dat, bias = b)))
})
print(round(bt, 3))
print(sapply(fits[c("SD", "SD_CE", "BHH_logCE")],
             function(f) mean(exp(residuals(f)))))

# Conventional thresholds for leverage and Cook's distance.
influential <- function(fit) {
  which(hatvalues(fit) > 2 * length(coef(fit)) / nobs(fit) |
          cooks.distance(fit) > 4 / nobs(fit))
}

# Identify high-leverage observations and those with large Cook's distance.
for (m in c("SD_CE", "SD_class", "BHH_logCE")) {
  f <- fits[[m]]
  infl <- data.frame(instance = dat$key, n = dat$n, type = dat$type,
                     hat = hatvalues(f), cook = cooks.distance(f))
  infl <- infl[influential(f), ]
  cat("\nInfluential instances:", m, "\n")
  print(infl[order(-infl$cook), ])
}
print(summary(MASS::rlm(models$SD_CE, data = dat)))
dfb <- dfbeta(fits$SD_CE)[, "log_sd:R"]
print(head(data.frame(instance = dat$key, n = dat$n, type = dat$type,
                      dfbeta = dfb)[order(-abs(dfb)), ], 10))

# Refit selected models after restricting size or excluding influential cases.
cl <- which(dat$type == "Clustered")
big <- head(cl[order(-dat$n[cl])], 5)
keep <- seq_len(nrow(dat))
cases <- list(All = keep, "n >= 50" = which(dat$n >= 50),
              "n >= 100" = which(dat$n >= 100), "n >= 200" = which(dat$n >= 200),
              "No 5 largest clustered" = setdiff(keep, big),
              "No influential (M2)" = setdiff(keep, influential(fits$SD_CE)),
              "No influential (M3)" = setdiff(keep, influential(fits$BHH_logCE)))
robustness <- lapply(names(cases), function(case) {
  d <- droplevels(dat[cases[[case]], ])
  cat("\nSensitivity:", case, "| n =", nrow(d), "\n")
  print(round(coef_table(lm(models$SD_CE, data = d)), 4))
  print(round(coef_table(lm(models$BHH_logCE, data = d)), 4))
  ed <- lapply(models[c("SD", "SD_CE", "BHH", "BHH_logCE")], repeated_cv, df = d)
  print(sapply(ed, function(e) sapply(groups(rowMeans(e), d$type), mean)))
  comparison <- comparison_table(ed, list(c("SD", "SD_CE"),
    c("BHH", "BHH_logCE"), c("SD_CE", "BHH_logCE")), d$type)
  print(comparison)
  data.frame(case = case, comparison)
})
write.csv(do.call(rbind, robustness), "results/tsplib_robustness.csv", row.names = FALSE)

# 5. Reliability of the Clark-Evans classification
# Compare the asymptotic test with Monte Carlo tests under CSR.
# Fewer simulations are used for very large instances to reduce computation.
set.seed(ANALYSIS_SEED)
dat$p_mc <- sapply(coords, function(cc) {
  pp <- as_ppp(cc)
  clarkevans.test(pp, correction = "Donnelly", method = "MonteCarlo",
                 nsim = if (pp$n <= 2000) 999 else 199)$p.value
})
dat$type_mc <- classify_ce(dat$R, dat$p_mc)
print(table(Asymptotic = dat$type, MonteCarlo = dat$type_mc))
write.csv(dat[, c("key", "n", "R", "p", "type", "p_mc", "type_mc")],
          "results/tsplib_classification.csv", row.names = FALSE)
dat$n_band <- cut(dat$n, c(0, 100, 1000, Inf),
                  labels = c("<=100", "101-1000", ">1000"))
print(table(dat$type, dat$n_band))
# Also classify by practical departure from CE=1, independent of the p-value.
for (delta in c(0.05, 0.10)) {
  lab <- ifelse(abs(dat$R - 1) < delta, "Random",
                ifelse(dat$R > 1, "Regular", "Clustered"))
  cat("\nEffect-size threshold:", delta, "\n")
  print(table(Test = dat$type, EffectSize = lab, dat$n_band))
}

# 6. Computing times
# Select 12 instances spanning the observed range of problem sizes.
# Compare the cost of CE, random-tour SD and distance-based benchmarks.
ix <- order(dat$n)[unique(round(seq(1, nrow(dat), length.out = 12)))]
cost <- do.call(rbind, lapply(ix, function(i) {
  cc <- coords[[i]]
  pp <- as_ppp(cc)
  # Measure CE 100 times per repetition and divide by 100.
  t_ce <- median(replicate(5, {
    system.time(for (j in 1:100) clarkevans(pp, correction = "Donnelly"))[["elapsed"]]
  })) / 100
  t_sd <- median(replicate(3, system.time(random_tour_sd(cc))[["elapsed"]]))
  t_nic <- system.time(nicola_features(cc))[["elapsed"]]
  data.frame(instance = dat$key[i], n = pp$n, t_CE_s = t_ce, t_SD_s = t_sd,
             t_Nicola_s = t_nic, CE_over_SD_pct = 100 * t_ce / t_sd)
}))
print(cost)
write.csv(cost, "results/computing_times.csv", row.names = FALSE)

# Instance-level instability of the Nicola benchmarks
# Inspect the upper tail of repeated-CV errors, not just average MAPE.
for (m in c("Nicola_back", "Nicola_fwd")) {
  inst <- data.frame(key = dat$key, n = dat$n, type = dat$type,
                     APE_mean = rowMeans(E[[m]]),
                     APE_median = apply(E[[m]], 1, median),
                     APE_max = apply(E[[m]], 1, max),
                     APE_M3 = rowMeans(E$BHH_logCE))
  print(head(inst[order(-inst$APE_mean), ], 10))
  print(sum(inst$APE_max > 50))
  print(rbind(Nicola = tapply(inst$APE_median, dat$type, mean),
              Model_3 = tapply(apply(E$BHH_logCE, 1, median), dat$type, mean)))
}

# Retain the analysis objects so script 06 can recreate the figures directly.
saveRDS(list(data = dat, errors = E, cv = tab_tsplib, paired = paired_tsplib,
             fits = fits, cost = cost, back_transform = bt), "results/tsplib.rds")

