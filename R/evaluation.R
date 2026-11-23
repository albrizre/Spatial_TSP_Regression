# Shared statistical functions for the simulation and TSPLIB analyses.
source("R/features.R")
# All models use the same CV and bootstrap settings.
CV_FOLDS <- 10L
CV_REPS <- 100L
BOOT_REPS <- 2000L
ANALYSIS_SEED <- 2026

# Cross-validation and paired inference.
# Absolute percentage error for each held-out instance.
ape <- function(pred, obs) 100 * abs(pred - obs) / obs

# Transform predictions back to tour-length units. Duan's smearing factor
# estimates E[exp(error)] from the training residuals only.
# Models fitted to the original response (rather than log(opt)) need no correction.
back_transform <- function(fit, newdata, bias = "smearing") {
  p <- predict(fit, newdata = newdata)
  if (identical(deparse(formula(fit)[[2]]), "opt")) return(p)
  switch(bias, naive = exp(p), lognormal = exp(p + summary(fit)$sigma^2 / 2),
         smearing = exp(p) * mean(exp(residuals(fit))))
}

# Repeated stratified cross-validation. Rows of E are instances; columns are
# repetitions. The seed is reset for every model, so all models receive the
# same sequence of folds. Model can be an lm formula or a fitting function.
repeated_cv <- function(model, df, k = CV_FOLDS, reps = CV_REPS,
                        bias = "smearing", seed = ANALYSIS_SEED) {
  set.seed(seed)
  E <- matrix(NA_real_, nrow(df), reps)

  for (r in seq_len(reps)) {
    # Allocate balanced random fold labels separately within each spatial class.
    folds <- integer(nrow(df))
    for (s in unique(df$type)) {
      i <- which(df$type == s)
      f <- rep(sample.int(k), length.out = length(i))
      folds[i] <- f[sample.int(length(i))]
    }

    # Fit on k-1 folds and predict the withheld observations.
    pred <- numeric(nrow(df))
    for (fold in unique(folds)) {
      train <- df[folds != fold, ]
      test <- folds == fold
      if (is.function(model)) {
        fit <- model(train)
      } else {
        fit <- lm(model, data = train)
      }
      # The smearing correction is estimated from training residuals only.
      pred[test] <- back_transform(fit, df[test, ], bias)
    }
    E[, r] <- ape(pred, df$opt)
  }
  E
}

# Report results separately by spatial class and for the full dataset.
groups <- function(x, type) c(split(x, type, drop = TRUE), list(Overall = x))

bootstrap_mean <- function(x, B) {
  replicate(B, mean(x[sample.int(length(x), replace = TRUE)]))
}

# Average the repeated APEs within each instance before computing MAPE.
# SD_reps measures variability across partitions; SE_boot resamples instances.
cv_table <- function(E, type, B = BOOT_REPS) {
  tab <- list()
  subsets <- groups(seq_along(type), type)
  for (m in names(E)) {
    inst <- rowMeans(E[[m]])
    for (g in names(subsets)) {
      i <- subsets[[g]]
      tab[[length(tab) + 1]] <- data.frame(
        model = m, type = g, n = length(i),
        MAPE = mean(inst[i]),
        SD_reps = sd(colMeans(E[[m]][i, , drop = FALSE])),
        SE_boot = sd(bootstrap_mean(inst[i], B))
      )
    }
  }
  do.call(rbind, tab)
}

# Paired comparisons use the average CV error of each instance.
# A positive difference favours the second model.
comparison_table <- function(E, pairs, type, B = BOOT_REPS) {
  tab <- list()
  for (pair in pairs) {
    d <- rowMeans(E[[pair[1]]]) - rowMeans(E[[pair[2]]])
    subsets <- groups(d, type)
    for (g in names(subsets)) {
      x <- subsets[[g]]
      ci <- quantile(bootstrap_mean(x, B), c(0.025, 0.975))
      tab[[length(tab) + 1]] <- data.frame(
        model1 = pair[1], model2 = pair[2], type = g, n = length(x),
        mean_diff = mean(x), median_diff = median(x),
        ci_low = unname(ci[1]), ci_high = unname(ci[2]),
        p_wilcoxon = if (length(x) >= 5) wilcox.test(x, exact = FALSE)$p.value else NA,
        share_B_better = mean(x > 0)
      )
    }
  }
  do.call(rbind, tab)
}

coef_table <- function(fit) {
  ct <- summary(fit)$coefficients
  ci <- confint(fit)[rownames(ct), , drop = FALSE]
  data.frame(Estimate = ct[, 1], SE = ct[, 2], CI_low = ci[, 1],
             CI_high = ci[, 2], p = ct[, 4])
}

fit_stats <- function(fit) {
  s <- summary(fit)
  c(adjR2 = s$adj.r.squared, sigma = s$sigma)
}

coefficient_table <- function(fits) {
  do.call(rbind, lapply(names(fits), function(m) {
    tab <- coef_table(fits[[m]])
    data.frame(model = m, term = rownames(tab), tab, row.names = NULL)
  }))
}
