# 06. Figures 1-4 from the manuscript (requires outputs from scripts 02 and 05).
# No models are refitted here; all plotted data come from the saved results.
library(ggplot2)
library(ggrepel)
dir.create("figures", showWarnings = FALSE)
tsplib <- readRDS("results/tsplib.rds")
simulation <- readRDS("results/simulation.rds")
dat <- tsplib$data
sim <- simulation$data
tab_tsplib <- tsplib$cv
tab_sim <- simulation$cv
pal <- c(Clustered = "#66c2a5", Random = "#fc8d62", Regular = "#8da0cb")
# Figure 1: Clark-Evans statistic by the simulated generating process.
p_ce <- ggplot(sim, aes(type, R, fill = type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1.5, position = position_jitter(width = 0.2, height = 0, seed = 2026)) +
  geom_hline(yintercept = 1, lty = "dashed", colour = "gray40") +
  scale_fill_manual(values = pal) + theme_bw() +
  labs(x = "Pattern type", y = "Clark-Evans statistic (CE)") +
  theme(legend.position = "none")
ggsave("figures/simulation_ce.pdf", p_ce, width = 7, height = 5)

# Figure 2: relationship between random-tour SD and exact optimal length.
p_sd <- ggplot(sim, aes(log_sd, log_opt, colour = type)) +
  geom_point() + scale_colour_manual(values = pal, name = "Pattern type") +
  theme_bw() + labs(x = "log(SD)", y = "log(L*)") +
  theme(legend.position = "bottom")
ggsave("figures/simulation_sd.pdf", p_sd, width = 7, height = 5)

# Figure 3: CE-based classification of the TSPLIB instances.
p_cls <- ggplot(dat, aes(n, R, colour = type)) +
  scale_x_log10() + geom_hline(yintercept = 1, lty = "dashed", col = "gray40") +
  geom_point() + theme_bw() +
  labs(title = "Classification of TSPLIB instances", x = "Number of nodes",
       y = "Clark-Evans statistic") +
  scale_colour_manual(values = pal, name = "Pattern type") +
  theme(text = element_text(size = 32), plot.title = element_text(face = "bold"),
        legend.position = "bottom")
p_cls <- p_cls + geom_label_repel(aes(label = key), box.padding = 0.35, seed = 2026,
                  point.padding = 0.5, segment.color = "grey50")
postscript("figures/tsplib_classification.eps", family = "Helvetica")
print(p_cls)
dev.off()

# Figure 4: out-of-sample MAPE for the main models on both datasets.
# Map corresponding simulation and TSPLIB models to common plot labels.
lab_tsplib <- c(SD = "Model 1", SD_CE = "Model 2", SD_class = "Model 1 by class",
                SD_classXCE = "Model 2 by class", BHH = "Scale only (BHH)",
                BHH_logCE = "Model 3")
lab_sim <- setNames(unname(lab_tsplib),
                    c("SD", "SD_CE", "SD_class", "SD_classXCE", "n_only", "n_logCE"))
plot_data <- function(tab, labels, dataset) {
  t <- tab[tab$model %in% names(labels), c("model", "type", "MAPE")]
  t$model <- unname(labels[t$model])
  data.frame(dataset = dataset, t)
}
oos <- rbind(plot_data(tab_sim, lab_sim, "Simulated"),
             plot_data(tab_tsplib, lab_tsplib, "TSPLIB"))
oos$model <- factor(oos$model, levels = unname(lab_tsplib))
oos$dataset <- factor(oos$dataset, levels = c("Simulated", "TSPLIB"))
oos$type <- factor(oos$type, levels = c("Clustered", "Random", "Regular", "Overall"))
write.csv(oos, "results/oos_mape_summary.csv", row.names = FALSE)
p_oos <- ggplot(oos, aes(model, MAPE, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~dataset, ncol = 1) + theme_bw() +
  scale_fill_manual(values = c(pal, Overall = "gray60"), name = "Pattern type") +
  labs(x = NULL, y = "Out-of-sample MAPE (%)") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("figures/oos_mape_summary.pdf", p_oos, width = 10, height = 7)
