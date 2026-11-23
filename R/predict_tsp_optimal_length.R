# Predict optimal Euclidean TSP tour lengths for new .tsp instances.
# Run source("R/predict_tsp_optimal_length.R") from the repository root.
source("R/features.R")

# instances: paths to Euclidean .tsp files containing NODE_COORD_SECTION.
# model: object returned by training_TSP_model(); if NULL, use the bundled
#        Model 3 fitted on 77 TSPLIB instances with known optimal lengths.
# conf_level: prediction-interval coverage (default 95%).
# Returns predictions in the original distance units of the supplied files.
predict_tsp_optimal_length <- function(instances, model = NULL, conf_level = 0.95,
                                       n_random = 1000,
                                       default_model_file = "models/default_tsp_model.rds") {
  if (is.null(model)) {
    if (!file.exists(default_model_file))
      stop("Default model not found: ", default_model_file,
           ". Run R from the repository root or supply a trained model.")
    model <- readRDS(default_model_file)
  }
  if (!inherits(model, "tsp_model") || !inherits(model$fit, "lm"))
    stop("'model' must be an object returned by training_TSP_model().")
  if (length(conf_level) != 1 || !is.numeric(conf_level) ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be a number between 0 and 1.")

  # Compute exactly the predictors required by the selected model.
  feats <- instance_features(instances, model$model, n_random)

  # predict.lm supplies t-based prediction intervals on the log scale.
  # Exponentiating gives intervals in tour-length units.
  pr <- predict(model$fit, newdata = feats, interval = "prediction",
                level = conf_level)
  out <- cbind(feats,
               fit = exp(pr[, "fit"]),
               fit_mean = exp(pr[, "fit"]) * model$smearing,
               lwr = exp(pr[, "lwr"]),
               upr = exp(pr[, "upr"]))

  # 'fit' is exp(predicted log length); 'fit_mean' additionally uses Duan's
  # smearing factor to correct retransformation bias in the estimated mean.
  # Flag extrapolation in the variables actually used by the model.
  r <- model$ranges
  outside <- rep(FALSE, nrow(out))
  if (model$model == "scale_logCE") {
    log_scale <- log(sqrt(out$n * out$A))
    outside <- log_scale < r$log_sqrt_nA[1] | log_scale > r$log_sqrt_nA[2] |
               out$CE < r$CE[1] | out$CE > r$CE[2]
  }
  if (model$model == "sd_ce") {
    outside <- out$CE < r$CE[1] | out$CE > r$CE[2]
  }
  if (tsp_model_needs_sd(model$model)) {
    if (is.null(r$sd)) stop("The supplied model has no training range for SD.")
    outside <- outside | out$sd < r$sd[1] | out$sd > r$sd[2]
  }
  out$extrapolation <- outside
  if (any(out$extrapolation))
    warning(sum(out$extrapolation), " instance(s) outside the training ranges.")

  rownames(out) <- NULL
  out
}
