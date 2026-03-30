library(TSP)
library(spatstat)
training_TSP_model <- function(training_instances) {
  
  ### Parameter description
  ### training_instances: A character vector containing the file paths to the .tsp files for training the regression models. The .tsp files should include the coordinates of the nodes under a NODE_COORD_SECTION section 

  df_analysis=c()
  for (i in 1:length(training_instances)){
    
    print(paste0("Instance number ",i))
    
    # Read instance i 
    aux <- read_TSPLIB(training_instances[i])
    if (!is.null(nrow(aux))){
      n=nrow(aux)
    } else{
      n=attr(aux,"Size")
    }
    
    # Solve instance by heuristic methods
    if (length(nrow(aux))>0){
      if (nrow(aux)>6000){
        tour <- solve_TSP(aux, method = "nn", two_opt = F, control = list(rep = 10))
      } else {
        tour <- solve_TSP(aux, method = "farthest_insertion", two_opt = F)
      }
    } else {
      if (length(aux)>6000){
        tour <- solve_TSP(aux, method = "nn", two_opt = F, control = list(rep = 10))
      } else {
        tour <- solve_TSP(aux, method = "farthest_insertion", two_opt = F)
      }
    }
    
    # Random tour generation for approximating the standard deviation of random tour length
    random_tours=c()
    for (k in 1:1000){
      tour_random <- solve_TSP(aux,method = "random")
      random_tours=c(random_tours,attr(tour_random,"tour_length"))
    }
    sd <- sd(random_tours)
    
    # Read coordinates
    lines <- readLines(training_instances[i])
    start_line <- grep("NODE_COORD_SECTION", lines) + 1
    end_line <- grep("EOF", lines) - 1
    if (length(end_line)==0){
      coords <- read.table(training_instances[i], skip = start_line - 1)
    } else {
      coords <- read.table(training_instances[i], skip = start_line - 1, nrows = end_line - start_line + 1)
    }
    if (ncol(coords)==3){coords=coords[,2:3]}
    
    # Clark-Evans statistic and classification of the pattern
    x_range <- range(coords[,1])
    y_range <- range(coords[,2])
    pp <- ppp(x = coords[,1],
              y = coords[,2],
              window = owin(x_range, y_range))
    result_test <- clarkevans.test(pp, correction = "Donnelly")
    clark_evans <- result_test$statistic
    if (result_test$p.value<0.05){
      if (clark_evans>1){
        type="regular"
      } else {
        type="clustered"
      }
    } else {
      type="random"
    }
    # Add data
    df_analysis=rbind(df_analysis,data.frame(n=n,opt=attr(tour,"tour_length"),sd=sd,CE=clark_evans,type=type))
  }
  rownames(df_analysis)=NULL
  
  model_coefficients=list()
  for (type in c("clustered","random","regular")){
    if (sum(df_analysis$type==type) >= 5){ # number of observations >= number of parameters to be estimated + 1
      lm_fit=lm(log(opt)~log(sd)*CE,data=df_analysis[df_analysis$type==type,])
      model_coefficients[[type]]=list(Intercept = as.numeric(lm_fit$coefficients[1]), 
                                      log_sd = as.numeric(lm_fit$coefficients[2]), 
                                      CE = as.numeric(lm_fit$coefficients[3]), 
                                      int = as.numeric(lm_fit$coefficients[4]), 
                                      rse = as.numeric(summary(lm_fit)$sigma))
    } else {
      model_coefficients[[type]]=list(Intercept = NA, 
                                      log_sd = NA, 
                                      CE = NA, 
                                      int = NA, 
                                      rse = NA)
    }
    
  }
  
  # Return coefficients
  return(model_coefficients)
}

