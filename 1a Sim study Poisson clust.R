library(TSP)
library(tspmeta)
library(spatstat)
library(eHDPrep)
library(spatstat)

# Set working directory

setwd("...")

# Instance simulation and covariates generation

N_random=1000
set.seed(12345)
data=c()
d_random_tours_list=list()
ecdf_random_tours_list=list()
count=0
data_list=list()
for (window_size in c(10)){
  for (n in c(10,20,50,100,200)){
    print(paste0("n = ",n))
    for (r in 1:10){
      print(r)
      count=count+1
      W=square(window_size)
      lambda=n/window_size^2
      cond=F
      while (cond==F){
        # sim_pattern=rCauchy(kappa=lambda/4,scale=0.02,mu=4,win=W) # kappa is the intensity of the parents, scale the cluster size, mu the number of offsprings
        sim_pattern=rCauchy(kappa=5/window_size^2,scale=0.25,mu=n/5,win=W,n.cond=n) # kappa is the intensity of the parents, scale the cluster size, mu the number of offsprings
        save(sim_pattern,file=paste0("Patterns/sim_pattern_clust_ins_",count,".rda"))
        result_test <- clarkevans.test(sim_pattern, correction = "Donnelly")
        if (sim_pattern$n>0){
          cond=T
        }
      }
      tsp=TSP(dist(data.frame(x=sim_pattern$x,y=sim_pattern$y)))
      # TSP solution
      TSP::write_TSPLIB(tsp,file=paste0("TSP_to_solve/clust_ins_",count,".tsp"))
      tour=solve_TSP(tsp)
      # Random tours
      n_points = length(sim_pattern$x)
      random_tours = numeric(N_random)
      tsp_mat <- as.matrix(tsp)
      for (i in 1:N_random){
        p <- sample(n_points)
        random_tours[i] <- sum(tsp_mat[cbind(p[-n_points], p[-1])]) + tsp_mat[p[n_points], p[1]]
      }
      data_list[[count]]=data.frame(count,window_size,n_exp=n,n_sim=length(sim_pattern$x),r,sd=sd(random_tours),opt=attr(tour,"tour_length"),
                                    R=result_test$statistic,p=result_test$p.value)
      # Empirical distribution tour length
      d_random_tours=density(random_tours)
      Af=approxfun(d_random_tours$x,d_random_tours$y)
      aux=Af(seq(0,2500,1))
      aux[is.na(aux)]=0
      d_random_tours_list[[count]]=cbind(x=seq(0,2500,1),y=aux)
      # ECDF
      ecdf_random_tours=ecdf(random_tours)
      aux=Af(seq(0,2500,1))
      aux=ecdf_random_tours(seq(0,2500,1))
      ecdf_random_tours_list[[count]]=cbind(x=seq(0,2500,1),y=aux)
    }
  }
}
data <- do.call(rbind, data_list)
data$Type="clust"
data_analysis=data[,c("sd","opt","Type","R")]
save(data,file="Results/data_clust.rda")
save(data_analysis,file="Results/data_analysis_clust.rda")
save(ecdf_random_tours_list,file="Results/ecdf_random_tours_list_clust.rda")
