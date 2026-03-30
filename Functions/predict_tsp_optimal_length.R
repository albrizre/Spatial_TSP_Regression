library(TSP)
library(spatstat)
predict_tsp_optimal_length <- function(instances, model_coefficients = NULL, conf_level = 0.95) {
  
  ### Parameter description
  ### instances: A character vector containing the file paths to the .tsp files for optimal tour length prediction. The .tsp files should include the coordinates of the nodes under a NODE_COORD_SECTION section 
  ### model_coefficients: A list including model coefficients for the pattern-specific regression models trained with function training_TSP_model
  ### conf_level: A value in ]0,1[ indicating the confidence level of the prediction interval to be obtained.
  
  df_results=c()
  for (i in 1:length(instances)){
    
    print(paste0("Instance number ",i))
    
    # Read instance i 
    aux <- read_TSPLIB(instances[i])
    if (!is.null(nrow(aux))){
      n=nrow(aux)
    } else{
      n=attr(aux,"Size")
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
    
  }
  
  for (i in 1:length(instances)){
    
    # Model coefficients
    if (is.null(training_instances)){
      params <- switch(tolower(type[i]),
                       "random"    = list(Intercept = 8.78539394, log_sd = -0.01831179, CE = -6.94703286, int = 0.94531265, rse = 0.1786966),
                       "clustered" = list(Intercept = -1.3943533, log_sd = 1.0934509, CE = 4.1000296, int = -0.1925425, rse = 0.2927208),
                       "regular"   = list(Intercept = 0.97385684,  log_sd = 0.98369281, CE = 0.51491472, int = -0.01522672, rse = 0.2327121),
                       stop("Invalid pattern type. Choose 'random', 'clustered', or 'regular'.")
      )
    } else {
      params <- switch(tolower(type[i]),
                       "random"    = list(Intercept = model_coefficients[["random"]]$Intercept, log_sd = model_coefficients[["random"]]$log_sd, 
                                          CE = model_coefficients[["random"]]$CE, int = model_coefficients[["random"]]$int, rse = model_coefficients[["random"]]$rse),
                       "clustered" = list(Intercept = model_coefficients[["clustered"]]$Intercept, log_sd = model_coefficients[["clustered"]]$log_sd, 
                                          CE = model_coefficients[["clustered"]]$CE, int = model_coefficients[["clustered"]]$int, rse = model_coefficients[["clustered"]]$rse),
                       "regular"   = list(Intercept = model_coefficients[["regular"]]$Intercept, log_sd = model_coefficients[["regular"]]$log_sd, 
                                          CE = model_coefficients[["regular"]]$CE, int = model_coefficients[["regular"]]$int, rse = model_coefficients[["regular"]]$rse),
                       stop("Invalid pattern type. Choose 'random', 'clustered', or 'regular'.")
      )
    }
    
    # Calculate Point Estimate (The Prediction)
    fit <- params$Intercept + (params$log_sd*log(sd[i])) + (params$CE*clark_evans[i]) + (params$int*log(sd[i])*clark_evans[i])
    
    # Calculate Prediction Interval
    # Note: This uses a simplified version assuming a large training sample (Z-score)
    # For exact T-distribution intervals, you'd need the degrees of freedom.
    alpha <- 1 - conf_level
    z_score <- qnorm(1 - alpha / 2)
    
    lower <- fit - (z_score * params$rse)
    upper <- fit + (z_score * params$rse)
    
    df_results=rbind(df_results,data.frame(
      instance = instances[i],
      sd = sd[i],
      clark_evans = clark_evans[i],
      type = type[i],
      fit = exp(fit),
      lwr = exp(lower),
      upr = exp(upper)
    ))
    
  }
  
  # Return results as a data.frame
  rownames(df_results) <- NULL
  return(df_results)
}


