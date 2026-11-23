# Fit a regression model for optimal Euclidean TSP tour lengths.
# Run source("R/training_TSP_model.R") from the repository root.
source("R/features.R")

# This helper is also used to fit Model 3 on an existing table of features.
# 'df' must contain the exact optimal length (opt), n, area (A) and CE;
# the SD-based models also need the standard deviation of random tours (sd).
make_tsp_model <- function(df, model = "scale_logCE") {
  model <- match.arg(model, names(tsp_model_formulas))
  required <- c("opt", "n", "A", "CE", if (tsp_model_needs_sd(model)) "sd")

  if (!all(required %in% names(df)))
    stop("Missing columns: ", paste(setdiff(required, names(df)), collapse = ", "))
  x <- as.matrix(df[, required, drop = FALSE])
  if (!is.numeric(x) || any(!is.finite(x)) || any(x <= 0))
    stop("Optimal lengths and required features must be positive and finite.")

  fit <- lm(tsp_model_formulas[[model]], data = df)
  if (anyNA(coef(fit)) || df.residual(fit) <= 0)
    stop("Not enough independent training information to estimate this model.")

  # Store the fitted regression, bias correction and training ranges.
  structure(list(
    model = model,
    fit = fit,
    smearing = mean(exp(residuals(fit))),
    n_train = nobs(fit),
    ranges = list(log_sqrt_nA = range(log(sqrt(df$n * df$A))),
                  CE = range(df$CE),
                  sd = if ("sd" %in% names(df)) range(df$sd)),
    train_data = df
  ), class = "tsp_model")
}

# training_instances: character vector with paths to Euclidean .tsp files.
# opt: exact optimal tour lengths in the same order as training_instances.
# model: "scale_logCE" (Model 3), "sd_ce" (Model 2), or "sd" (Model 1).
# Returns an object that can be passed to predict_tsp_optimal_length().
training_TSP_model <- function(training_instances, opt, model = "scale_logCE",
                               n_random = 1000) {
  model <- match.arg(model, names(tsp_model_formulas))
  if (!is.character(training_instances) || !length(training_instances))
    stop("Provide paths to the training .tsp files.")
  if (missing(opt) || !is.numeric(opt) || length(opt) != length(training_instances) ||
      any(!is.finite(opt)) || any(opt <= 0))
    stop("Provide one positive, finite optimal tour length per training instance.")

  # Compute features from the .tsp coordinates; Model 3 skips random tours.
  df <- instance_features(training_instances, model, n_random)
  df$opt <- opt
  make_tsp_model(df, model)
}

print.tsp_model <- function(x, ...) {
  cat("TSP model:", x$model, "| Training instances:", x$n_train, "\n")
  print(summary(x$fit)$coefficients)
  cat("Smearing factor:", round(x$smearing, 4), "\n")
  invisible(x)
}
