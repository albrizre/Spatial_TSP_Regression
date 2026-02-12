library(TSP)
library(tspmeta)

# Set working directory

setwd("...")

# Read instances from TSPlib. For reproducibility, the user should have a Data/ALL_tsp/ folder with all TSPLIB instances in .tsp format
# These instances can be obtained from https://github.com/mastqe/tsplib

instances=dir("Data/ALL_tsp/")
instances=instances[-grep("tour",instances)]

set.seed(123)
data=c()
ecdf_random_tours_list=c()
# This process takes several hours given the size of some of the instances
for (instance in instances){

  print(instance)

  skip_to_next <- FALSE
  tryCatch(aux <- read_TSPLIB(paste0("Data/ALL_tsp/",instance,"/",instance)), error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) {
    print("Error")
    next
    } else {
    print(Sys.time())
    if (length(nrow(aux))>0){
      if (nrow(aux)>10000){
        tour <- solve_TSP(aux, method = "nn", two_opt = F, control = list(rep = 10))
      } else {
        tour <- solve_TSP(aux, method = "farthest_insertion", two_opt = F)
      }
    } else {
      if (length(aux)>10000){
        tour <- solve_TSP(aux, method = "nn", two_opt = F, control = list(rep = 10))
      } else {
        tour <- solve_TSP(aux, method = "farthest_insertion", two_opt = F)
      }
    }
    print(Sys.time())
    if (!is.null(nrow(aux))){
      n=nrow(aux)
    } else{
      n=attr(aux,"Size")
    }
    random_tours=c()
    for (i in 1:1000){
      print(i)
      tour_random <- solve_TSP(aux,method = "random")
      random_tours=rbind(random_tours,data.frame(id=i,Length=attr(tour_random,"tour_length")))
    }
    data=rbind(data,data.frame(instance,n=n,sd=sd(random_tours$Length),opt=attr(tour,"tour_length")))
    # Empirical distribution tour length
    d_random_tours=density((random_tours$Length-min(random_tours$Length))/(max(random_tours$Length)-min(random_tours$Length)))
    Af=approxfun(d_random_tours$x,d_random_tours$y)
    aux=Af(seq(0,1,0.01))
    aux[is.na(aux)]=0
    # ECDF
    ecdf_random_tours=ecdf((random_tours$Length-min(random_tours$Length))/(max(random_tours$Length)-min(random_tours$Length)))
    aux=Af(seq(0,1,0.01))
    aux=ecdf_random_tours(seq(0,1,0.01))
    ecdf_random_tours_list[[instance]]=cbind(x=seq(0,1,0.01),y=aux)
  }
}
save(data,file="Results/TSPlib_data_farthest_insertion.rda")
save(ecdf_random_tours_list,file="Results/TSPlib_ecdfs.rda")

plot(log(data$sd),log(data$opt))
cor(log(data$sd),log(data$opt))


