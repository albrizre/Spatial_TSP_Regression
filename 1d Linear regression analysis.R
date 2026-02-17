library(ggplot2)
library(DAAG)

# Set working directory

setwd("...")

# Dataset construction

# We solved the simulated instances with Concorde externally. Thus, we provide an Excel file with the optimal tour lengths obtained as an extra file
# Another option would be solving them with TSP_solve, calling the Concorde solver. In this case, generate_final_data.R would not be required,
# and optimal tour lengths could be obtained by adapting codes 1a, 1b, and 1c

source("Functions/generate_final_data.R")
data=data[-which(data$sd<1),]

df_plot=data
df_plot$Type[df_plot$Type=="clust"]="Clustered"
df_plot$Type[df_plot$Type=="reg"]="Regular"
df_plot$Type[df_plot$Type=="hom"]="Random"
ggplot(df_plot, aes(x = Type, y = R, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 1, size = 1.5) +
  theme_bw() +
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"),name="Pattern type")+
  geom_hline(yintercept = 1,lty="dashed",col="gray40")+
  labs(
    title = "Distribution of CE by pattern type",
    x = "Pattern type",
    y = "Clark-Evans statistic (CE)",
    fill = "Type"
  ) +
  theme(text = element_text(size=32),
        legend.position = "none",
        plot.title = element_text(face="bold"))

# Global model fit

lm_fit_global=lm(log(opt)~log(sd),data=data)
summary(lm_fit_global)
lm_fit_global_R=lm(log(opt)~log(sd)*R,data=data)
summary(lm_fit_global_R)
print(mean(100*abs(exp(lm_fit_global$fitted.values[data$Type=="clust"])-exp(lm_fit_global$model[data$Type=="clust",1]))/exp(lm_fit_global$model[data$Type=="clust",1])))
print(mean(100*abs(exp(lm_fit_global_R$fitted.values[data$Type=="clust"])-exp(lm_fit_global_R$model[data$Type=="clust",1]))/exp(lm_fit_global_R$model[data$Type=="clust",1])))
print(mean(100*abs(exp(lm_fit_global$fitted.values[data$Type=="hom"])-exp(lm_fit_global$model[data$Type=="hom",1]))/exp(lm_fit_global$model[data$Type=="hom",1])))
print(mean(100*abs(exp(lm_fit_global_R$fitted.values[data$Type=="hom"])-exp(lm_fit_global_R$model[data$Type=="hom",1]))/exp(lm_fit_global_R$model[data$Type=="hom",1])))
print(mean(100*abs(exp(lm_fit_global$fitted.values[data$Type=="inhib"])-exp(lm_fit_global$model[data$Type=="inhib",1]))/exp(lm_fit_global$model[data$Type=="inhib",1])))
print(mean(100*abs(exp(lm_fit_global_R$fitted.values[data$Type=="inhib"])-exp(lm_fit_global_R$model[data$Type=="inhib",1]))/exp(lm_fit_global_R$model[data$Type=="inhib",1])))
save(lm_fit_global,file="Results/lm_fit_global.rda")
save(lm_fit_global_R,file="Results/lm_fit_global_R.rda")

# Interaction analysis (ANCOVA) and scatter plot

lm_fit_int=lm(log(data$opt)~log(data$sd)*data$Type)
summary(lm_fit_int)

df_plot=data
df_plot$Type[df_plot$Type=="clust"]="Clustered"
df_plot$Type[df_plot$Type=="reg"]="Regular"
df_plot$Type[df_plot$Type=="hom"]="Random"
ggplot(data=df_plot,aes(x=log(sd),y=log(opt),col=Type))+
  geom_point(size=2)+
  theme_bw() +
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb"),name="Pattern type")+
  labs(
    title = paste("Scatter plot by pattern type"),
    y = bquote(log(italic(L)^"*")),
    x = expression(log(italic(SD)))
  ) +
  theme(text = element_text(size=32),
        legend.position = "bottom",
        plot.title = element_text(face="bold"))

# Pattern-specific model fit

lm_fit_hom=lm(log(opt)~log(sd),data=data,subset = (Type == "hom"))
lm_fit_clust=lm(log(opt)~log(sd),data=data,subset = (Type == "clust"))
lm_fit_reg=lm(log(opt)~log(sd),data=data,subset = (Type == "reg"))
lm_fit_hom_R=lm(log(opt)~log(sd)*R,data=data,subset = (Type == "hom"))
lm_fit_clust_R=lm(log(opt)~log(sd)*R,data=data,subset = (Type == "clust"))
lm_fit_reg_R=lm(log(opt)~log(sd)*R,data=data,subset = (Type == "reg"))
summary(lm_fit_clust)
summary(lm_fit_hom)
summary(lm_fit_reg)
round(confint(lm_fit_clust),2)
round(confint(lm_fit_hom),2)
round(confint(lm_fit_reg),2)
summary(lm_fit_clust_R)
summary(lm_fit_hom_R)
summary(lm_fit_reg_R)
save(lm_fit_hom,file="Results/lm_fit_hom.rda")
save(lm_fit_clust,file="Results/lm_fit_clust.rda")
save(lm_fit_reg,file="Results/lm_fit_reg.rda")
save(lm_fit_hom_R,file="Results/lm_fit_hom_R.rda")
save(lm_fit_clust_R,file="Results/lm_fit_clust_R.rda")
save(lm_fit_reg_R,file="Results/lm_fit_reg_R.rda")

# MAPE for pattern-specific models (without and with Clark-Evans)
for (type in c("clust","hom","reg")){

  aux=data[data$Type==type,]
  aux$log_sd=log(aux$sd)
  aux$log_opt=log(aux$opt)

  lm_fit<-lm(log_opt~log_sd,data=aux)
  print(mean(100*abs(exp(lm_fit$fitted.values)-exp(lm_fit$model[,1]))/exp(lm_fit$model[,1])))

  lm_fit<-lm(log_opt~log_sd*R,data=aux)
  print(mean(100*abs(exp(lm_fit$fitted.values)-exp(lm_fit$model[,1]))/exp(lm_fit$model[,1])))

}




