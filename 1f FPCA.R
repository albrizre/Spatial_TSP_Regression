library(TSP)
library(tspmeta)
library(spatstat)
library(eHDPrep)
library(fda)
library(GGally)
library(fdapace)
library(tidyr)

# Set working directory

setwd("...")

# Loading ecdfs

load("Results/ecdf_random_tours_list_hom.rda")
ecdf_random_tours_list_all=ecdf_random_tours_list[1:50]
load("Results/ecdf_random_tours_list_clust.rda")
ecdf_random_tours_list_all=c(ecdf_random_tours_list_all,ecdf_random_tours_list[1:50])
load("Results/ecdf_random_tours_list_reg.rda")
ecdf_random_tours_list_all=c(ecdf_random_tours_list_all,ecdf_random_tours_list[1:50])

# FPCA

argvals <- ecdf_random_tours_list_all[[1]][, 1]
data_matrix <- sapply(ecdf_random_tours_list_all, function(x) x[, 2])
basis <- create.bspline.basis(rangeval = range(argvals), nbasis = 20)
ecdf_fd <- smooth.basis(argvals, data_matrix, basis)$fd
fpca_results <- pca.fd(ecdf_fd, nharm = 5)
# plot.pca.fd(fpca_results, ask = FALSE)

data_fpca=data.frame(SC1=fpca_results$scores[,1],SC2=fpca_results$scores[,2],
                     SC3=fpca_results$scores[,3],SC4=fpca_results$scores[,4],
                     Type=c(rep("Random",50),rep("Clustered",50),rep("Regular",50)),
                     n=data$n_sim)

print(ggplot(data_fpca, aes(x = SC1, y = SC2, color = Type)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "FPC1 vs. FPC2 scores")+
  xlab("FPC1 score")+ylab("FPC2 score")+
  theme(text = element_text(size=32),
        plot.title = element_text(face="bold"),
        legend.position = "bottom")+
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb"),name="Pattern type"))

# Plot main eigenfunctions

eval_grid <- seq(min(argvals), max(argvals), length.out = 500)
harm_eval <- eval.fd(eval_grid, fpca_results$harmonics[1:3])
harm_df <- data.frame(
  x = eval_grid,
  FPC1 = harm_eval[, 1],
  FPC2 = harm_eval[, 2]
) %>%
  pivot_longer(cols = c(FPC1, FPC2), names_to = "Eigenfunction", values_to = "Value")
ggplot(harm_df, aes(x = x, y = Value, color = Eigenfunction)) +
  geom_line(linewidth = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 1) +
  labs(
    title = "FPC Eigenfunctions",
    x = expression(paste("Random tour length (", italic(L), ")")),
    y = "Value"
  ) +
  theme_bw() +
  theme(text = element_text(size=32),
        plot.title = element_text(face="bold"),
        legend.position = "bottom")+
  scale_color_manual(values = c("FPC1" = "#01665e", "FPC2" = "#80cdc1")) +
  theme(legend.position = "bottom")

