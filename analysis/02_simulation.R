# 02. Analyse the main simulated dataset (Sections 3.1.1-3.1.2).
# Input: archived features and exact lengths; output: tables and fitted models.
source("R/setup.R")
source("R/simulation_data.R")
sim <- load_simulation()
print(table(sim$type))

# Compare Models 1 and 2, class-specific versions and alternative scale
# predictors. Since A is fixed, Model 3 reduces to log(n) + log(CE).
models_sim <- list(
  SD = log_opt ~ log_sd,
  SD_CE = log_opt ~ log_sd * R,
  SD_class = log_opt ~ log_sd * type,
  SD_class_CE = log_opt ~ log_sd * type + log_sd * R,
  SD_classXCE = log_opt ~ log_sd * R * type,
  n_only = log_opt ~ log_n,
  n_CE = log_opt ~ log_n + R,
  n_logCE = log_opt ~ log_n + log(R),
  SD_n_CE = log_opt ~ log_sd + log_n + R
)
# Out-of-sample errors: the same stratified folds are used for every model.
E_sim <- list()
for (m in names(models_sim)) {
  E_sim[[m]] <- repeated_cv(models_sim[[m]], sim)
}
tab_sim <- cv_table(E_sim, sim$type)
paired_sim <- comparison_table(E_sim,
  list(c("SD", "SD_CE"), c("SD_class", "SD_class_CE"),
       c("SD_CE", "SD_class_CE"), c("n_only", "n_logCE"),
       c("SD", "SD_n_CE"), c("n_logCE", "SD_n_CE")), sim$type)
print(tab_sim)
print(paired_sim)
write.csv(tab_sim, "results/simulation_cv.csv", row.names = FALSE)
write.csv(paired_sim, "results/simulation_paired.csv", row.names = FALSE)

# Full-sample fits for inference on regression coefficients, not CV accuracy.
fits_sim <- lapply(models_sim, lm, data = sim)
print(sapply(fits_sim, fit_stats))
print(coef_table(fits_sim$SD_CE))
print(coef_table(fits_sim$n_logCE))
write.csv(coefficient_table(fits_sim), "results/simulation_coefficients.csv", row.names = FALSE)
# Fit the SD baseline separately within each generating pattern class.
class_fits <- lapply(levels(sim$type), function(type) {
  cat("\nModel 1:", type, "\n")
  f <- lm(log_opt ~ log_sd, data = sim[sim$type == type, ])
  tab <- coef_table(f)
  print(cbind(tab, t = summary(f)$coefficients[, 3]))
  data.frame(type = type, term = rownames(tab), tab,
             t = summary(f)$coefficients[, 3], row.names = NULL)
})
write.csv(do.call(rbind, class_fits), "results/simulation_class_coefficients.csv",
          row.names = FALSE)
# Percentage of correct CE-based classifications by target n and true type.
classification <- 100 * tapply(classify_ce(sim$R, sim$p) == sim$type,
                               list(sim$n_exp, sim$type), mean)
print(round(classification, 1))
write.csv(classification, "results/simulation_classification.csv")
# Save objects for later figures and additional analyses.
saveRDS(list(data = sim, errors = E_sim, cv = tab_sim, paired = paired_sim,
             fits = fits_sim), "results/simulation.rds")

