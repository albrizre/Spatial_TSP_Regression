library(ggplot2)
library(DAAG)
library(fda.usc)
library(fda)

# Set working directory

setwd("...")

df_plot=c()
df_error=c()
results_table=c()
for (type in c("hom","reg","clust")){

  load(paste0("Results/lm_fit_",type,".rda"))
  if (type=="hom"){subtitle_aux="Homogeneous instance";lm_fit=lm_fit_hom}
  if (type=="reg"){subtitle_aux="Regular instance";lm_fit=lm_fit_reg}
  if (type=="clust"){subtitle_aux="Clustered instance";lm_fit=lm_fit_clust}

  load(paste0("Results/ecdf_random_tours_list_",type,".rda"))

  source("Functions/generate_final_data.R")
  data=data[-which(data$sd<1),]
  data=data[data$Type==type,]

  # load(paste0("Results/data_",type,".rda"))

  # FDA regression: ecdf

  matrix_ecdf=c()
  matrix_x1=c()
  for (i in which(data$sd>=1)){
    matrix_x1=rbind(matrix_x1,ecdf_random_tours_list[[i]][,1])
    matrix_ecdf=rbind(matrix_ecdf,ecdf_random_tours_list[[i]][,2])
  }
  dim(matrix_ecdf)
  dim(matrix_x1)

  data_fda = fdata(matrix_ecdf, argvals = matrix_x1[1,])
  rangett = c(0,2500)
  k1=40
  k2=4
  basis1 = create.bspline.basis(rangeval = rangett, nbasis = k1) # for the covariate (we can use large k)
  basis2 = create.bspline.basis(rangeval = rangett, nbasis = k2) # for the covariate-response curve (better a low k for regularization)
  reg = fregre.basis(data_fda, log(data$opt[which(data$sd>=1)]), basis.x = basis1, basis.b = basis2)

  matrix_interaction = matrix_ecdf*data$R
  data_fda_int = fdata(matrix_interaction, argvals = matrix_x1[1,])
  datlist = list(
    df = data.frame(y = log(data$opt[which(data$sd>=1)]), R = data$R),
    X = data_fda,
    X_int = data_fda_int
  )
  reg_int = fregre.lm(y ~ R + X + X_int, data = datlist,
                  basis.x = list(X = basis1, X_int = basis1),
                  basis.b = list(X = basis2, X_int = basis2))
  results_table=rbind(results_table,data.frame(type,
                                               MAPE=mean(100*abs(reg$fitted.values-reg$y)/reg$y),
                                               MAPE_int=mean(100*abs(reg_int$fitted.values-log(data$opt[which(data$sd>=1)]))/log(data$opt[which(data$sd>=1)])),
                                               R2=reg$r2,
                                               R2_int=reg_int$r2))

  # Bootstrap for confidence band

  argvals_grid <- data_fda$argvals
  B <- 100
  beta_samples <- matrix(NA, nrow = B, ncol = length(argvals_grid))
  set.seed(123)
  for (b in 1:B) {
    # print(b)
    idx <- sample(1:nrow(matrix_ecdf), replace = TRUE)
    boot_fda <- data_fda[idx,]
    boot_y <- log(data$opt[idx])
    boot_reg <- fregre.basis(boot_fda, boot_y, basis.x = basis1, basis.b = basis2)
    beta_samples[b,] <- as.numeric(eval.fd(argvals_grid, boot_reg$beta.est))
  }
  beta_mean <- colMeans(beta_samples)
  beta_lower <- apply(beta_samples, 2, quantile, probs = 0.025)
  beta_upper <- apply(beta_samples, 2, quantile, probs = 0.975)

  df_plot <- rbind(df_plot,data.frame(type, x = argvals_grid, beta = beta_mean, lower = beta_lower, upper = beta_upper))

  df_comp <- data.frame(
    Observed = rep(log(data$opt[which(data$sd>=1)]), 2),
    Fitted   = c(lm_fit$fitted.values, reg$fitted.values),
    Model    = rep(c("Linear regression model (SD)", "Functional regression model (CDF)"), each = length(lm_fit$fitted.values))
  )

  # Scatter plot observed vs. fitted
  print(ggplot(df_comp, aes(x = Observed, y = Fitted, color = Model)) +
    geom_point() + # Alpha helps see overlapping points
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    # geom_smooth(method = "lm", se = FALSE, linetype = "solid") + # Adds a trend line for each
    labs(
      title = "Observed vs. predicted optimal values",
      subtitle = subtitle_aux,
      x = "Observed optimal value",
      y = "Predicted optimal value"
    ) +
    scale_color_manual(values = c("#6baed6","#08519c"))+
    theme_bw()+
    theme(text = element_text(size = 24),
    axis.title = element_text(face = "bold"),
    title = element_text(colour = "black",face = "bold"),
    legend.position = "bottom"))

  obs = log(data$opt)
  error_lm  <- 100*abs(obs - lm_fit$fitted.values)/obs
  error_fda <- 100*abs(obs - reg$fitted.values)/obs
  df_error <- rbind(df_error,data.frame(type,
    Error = c(error_lm, error_fda),
    Model = rep(c("Linear regression model (SD)", "Functional regression model (CDF)"), each = length(obs))))


}

# Plot beta(L) effect
df_plot=df_plot[df_plot$x<=1200,]
df_plot$type[df_plot$type=="clust"]="Clustered"
df_plot$type[df_plot$type=="reg"]="Regular"
df_plot$type[df_plot$type=="hom"]="Random"
line_colors <- c("#66c2a5", "#fc8d62", "#8da0cb")
fill_colors <- c("#d1ede4", "#fed8c8", "#dde2f0")
ggplot(df_plot, aes(x = x, col = type, fill = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA) +
  geom_line(aes(y = beta), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~type) +
  theme_bw() +
  scale_color_manual(values = line_colors, name = "Pattern type") +
  scale_fill_manual(values = fill_colors, name = "Pattern type") +
  scale_x_continuous(breaks=seq(0,800,400))+
  labs(
    title = "Estimated length-dependent covariate effects",
    y = expression(beta(italic(L))),
    x = expression(paste("Random tour length (", italic(L), ")"))) +
  theme(
    text = element_text(size = 24),
    axis.title = element_text(face = "bold"),
    title = element_text(colour = "black", face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none")

# Plot error distribution
df_error$Model=factor(df_error$Model,levels=c("Linear regression model (SD)", "Functional regression model (CDF)"))
df_error$type[df_error$type=="clust"]="Clustered"
df_error$type[df_error$type=="reg"]="Regular"
df_error$type[df_error$type=="hom"]="Random"
print(ggplot(df_error, aes(x = log(Error), col = Model)) +
        geom_density(size=1) + # Alpha makes the overlap visible
        facet_wrap(~type) +
        theme_bw() +
        labs(
          title = "Error distribution comparison",
          x = "Percent absolute error (log-transformed)",
          y = "Density"
        ) +
        scale_color_manual(values = c("#6baed6", "#08519c")) +
        theme_bw()+
        theme(text = element_text(size = 24),
              axis.title = element_text(face = "bold"),
              title = element_text(colour = "black",face = "bold"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              legend.position = "bottom"))

# View results
results_table[c(3,1,2),]
