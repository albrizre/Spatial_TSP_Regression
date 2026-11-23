# Common setup for the numbered analysis scripts.
# Fix the random-number generator and seed for reproducible model comparisons.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(2026)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
source("R/evaluation.R")
# Record the software environment used for the analysis.
capture.output(sessionInfo(), file = "results/session_info.txt")
