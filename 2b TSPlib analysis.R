library(ggplot2)
library(TSP)
library(tspmeta)
library(spatstat)
library(ggrepel)

# Set working directory

setwd("...")

# Data loading

load("Results/TSPlib_data_farthest_insertion.rda")

# Classification of instances

instances=data$instance
data_class=c()
for (instance in instances){

  print(instance)

  coord_aux=openxlsx::read.xlsx("TSPlib_coords.xlsx",sheet = instance, colNames = F)

  if (!is.null(nrow(coord_aux))){
    coords=c()
    for (i in 1:nrow(coord_aux)){
      print(i)
      split_row=unlist(strsplit(coord_aux[i,]," "))
      split_row=split_row[split_row!=""]
      if (length(split_row)>2){split_row=split_row[2:3]}
      coords=rbind(coords,as.numeric(split_row))
    }

    x_range <- range(coords[,1])
    y_range <- range(coords[,2])
    pp <- ppp(x = coords[,1],
              y = coords[,2],
              window = owin(x_range, y_range))
    result_test <- clarkevans.test(pp, correction = "Donnelly")

    ancho_original <- diff(range(pp$x))
    alto_original <- diff(range(pp$y))
    max_lado <- max(ancho_original, alto_original)
    pp_norm <- rescale(pp, max_lado / 10)

    if (result_test$p.value<0.05){
      if (result_test$statistic>1){
        type="regular"
        mu=NA;kappa=NA;scale=NA;mu_norm=NA;kappa_norm=NA;scale_norm=NA
      } else {
        type="clustered"
        fit_cauchy <- kppm(pp ~ 1, clusters = "Cauchy")
        mu=fit_cauchy$mu;kappa=fit_cauchy$clustpar[1];scale=fit_cauchy$clustpar[2]
        fit_cauchy_norm <- kppm(pp_norm ~ 1, clusters = "Cauchy")
        mu_norm=fit_cauchy_norm$mu;kappa_norm=fit_cauchy_norm$clustpar[1];scale_norm=fit_cauchy_norm$clustpar[2]
      }
    } else {
      type="random"
      mu=NA;kappa=NA;scale=NA;mu_norm=NA;kappa_norm=NA;scale_norm=NA
    }

    data_class=rbind(data_class,data.frame(instance,n=nrow(coord_aux),R=result_test$statistic,p=result_test$p.value,type=type,
                                           mu,kappa,scale,mu_norm,kappa_norm,scale_norm))
  }

}
save(data_class,file="Results/TSPLib_data_class.rda")

load(file="Results/TSPLib_data_class.rda")
data_class$instance2=gsub(".tsp","",data_class$instance)
data_class$type[data_class$type=="clustered"]="Clustered"
data_class$type[data_class$type=="random"]="Random"
data_class$type[data_class$type=="regular"]="Regular"
ggplot(data_class, aes(x = n, y = R, col = type)) +
  # Transform x axis to log scale
  scale_x_log10() +
  # Use geom_label_repel to show instance names without overlapping
  geom_label_repel(aes(label = instance2),
                   box.padding = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  # Adding a point as well so the label has a clear "anchor"
  geom_point() +
  theme_bw() +
  labs(
    title = "Classification of TSPLIB instances",
    x = "Number of nodes",
    y = "Clark-Evans statistic",
  )+
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb"),name="Pattern type")+
  theme(text = element_text(size=32),
        plot.title = element_text(face="bold"),
        legend.position = "bottom")

# Data analyis

data_analysis=merge(data,data_class,by="instance")

lm_fit_global=lm(log(opt)~log(sd),data=data_analysis)
lm_fit_global_R=lm(log(opt)~log(sd)*R,data=data_analysis)

print(round(summary(lm_fit_global_R)$adj.r.squared,2))
print(round(mean(100*abs(exp(lm_fit_global$fitted.values[data_analysis$type=="Clustered"])-exp(lm_fit_global$model[data_analysis$type=="Clustered",1]))/exp(lm_fit_global$model[data_analysis$type=="Clustered",1])),2))
print(round(mean(100*abs(exp(lm_fit_global$fitted.values[data_analysis$type=="Random"])-exp(lm_fit_global$model[data_analysis$type=="Random",1]))/exp(lm_fit_global$model[data_analysis$type=="Random",1])),2))
print(round(mean(100*abs(exp(lm_fit_global$fitted.values[data_analysis$type=="Regular"])-exp(lm_fit_global$model[data_analysis$type=="Regular",1]))/exp(lm_fit_global$model[data_analysis$type=="Regular",1])),2))

print(round(summary(lm_fit_global_R)$adj.r.squared,2))
print(round(mean(100*abs(exp(lm_fit_global_R$fitted.values[data_analysis$type=="Clustered"])-exp(lm_fit_global_R$model[data_analysis$type=="Clustered",1]))/exp(lm_fit_global_R$model[data_analysis$type=="Clustered",1])),2))
print(round(mean(100*abs(exp(lm_fit_global_R$fitted.values[data_analysis$type=="Random"])-exp(lm_fit_global_R$model[data_analysis$type=="Random",1]))/exp(lm_fit_global_R$model[data_analysis$type=="Random",1])),2))
print(round(mean(100*abs(exp(lm_fit_global_R$fitted.values[data_analysis$type=="Regular"])-exp(lm_fit_global_R$model[data_analysis$type=="Regular",1]))/exp(lm_fit_global_R$model[data_analysis$type=="Regular",1])),2))

for (type in c("Clustered","Random","Regular")){

  aux=data_analysis[data_analysis$type==type,]
  lm_fit<-lm(log(opt)~log(sd),data=aux)
  lm_fit_R<-lm(log(opt)~log(sd)*R,data=aux)
  print(round(mean(100*abs(exp(lm_fit$fitted.values)-exp(lm_fit$model[,1]))/exp(lm_fit$model[,1])),2))
  print(round(mean(100*abs(exp(lm_fit_R$fitted.values)-exp(lm_fit_R$model[,1]))/exp(lm_fit_R$model[,1])),2))
  print(round(summary(lm_fit)$adj.r.squared,2))
  print(round(summary(lm_fit_R)$adj.r.squared,2))

}

# Table with the classification

data_class$p=round(data_class$p,2)
data_class$R=round(data_class$R,2)
data_class$instance=gsub(".tsp","",data_class$instance)
df_subset <- data_class[, c("instance", "n", "R", "p", "type")]
n_rows <- nrow(df_subset)
half <- ceiling(n_rows / 2)
part1 <- df_subset[1:half, ]
part2 <- df_subset[(half + 1):n_rows, ]
if (nrow(part2) < nrow(part1)) {
  part2[nrow(part1), ] <- ""
}
combined_df <- cbind(part1, part2)
library(xtable)
print(xtable(combined_df), include.rownames = FALSE)


