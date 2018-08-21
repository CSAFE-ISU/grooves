cp1_gibbs <- function(data, iter, start.vals, prop_var, cp_prop_var, tol = 10, warmup = 5000, verbose = FALSE)
{
  ##data is a data frame with column x and column y
  
  ## initialize parameter list
  par <- list()
  par$sigma <- matrix(nrow = warmup + 1, ncol = 2) ## the variance of the GP
  par$sigma[1,] <- start.vals$sigma
  
  par$l <- matrix(nrow = warmup + 1, ncol = 2) ## length scale of the GP
  par$l[1,] <- start.vals$l
  
  par$tau <- matrix(nrow = warmup + 1, ncol = 2) ## nugget of the data model
  par$tau[1,] <- start.vals$tau
  
  par$cp <- matrix(nrow = warmup + 1, ncol = 1) ## changepoint locations
  par$cp[1,] <- start.vals$cp
  
  ## range on the x-axis of data
  interval <- range(data$x)
  
  ## current values of parameters
  sigma <- start.vals$sigma
  l <- start.vals$l
  tau <- start.vals$tau
  cp <- start.vals$cp
  
  ## initialize acceptance rates
  accept <- list()
  accept$gp_par <- matrix(data = c(0,0), nrow = 1, ncol = 2)
  accept$cp <- 0
  
  ## gibbs warmup iterations
  for(i in 1:(warmup))
  {
    xrange <- matrix(nrow = 2, ncol = 2)
    xrange[1,] <- c(interval[1], cp[1])
    xrange[2,] <- c(cp[1], interval[2])
    
    ## given changepoints make proposal for MH steps for GP parameters
    for(j in 1:2)
    {
      prop <- as.numeric(rmvnorm(n = 1, mean = c(sigma[j], l[j], tau[j]), sigma = prop_var[[j]]))
      if(verbose == TRUE)
      {
        print(paste("iteration: ",i))
        print(paste(j,"-th GP parameter proposal: ", prop))
      }
      
      ## skip this chunk of data if the proposal produces negative values
      if(any(prop <= 0))
      {
        # par$sigma[i,j] <- sigma[j]
        # par$l[i,j] <- l[j]
        # par$tau[i,j] <- tau[j]
        next
      }
      else{
        temp_dat <- data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$y
        
        ## proposal doesn't appear because it should cancel
        log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = prop[1], l = prop[2]) +
          dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) +
          dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) -
          (lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma[j], l = l[j]) +
             dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) +
             dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE))
        
        if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
        {
          sigma[j] <- prop[1]
          l[j] <- prop[2]
          tau[j] <- prop[3]
        }
      }
    }
    ## update GP parameters
    par$sigma[i + 1,] <- sigma
    par$l[i + 1,] <- l
    par$tau[i + 1,] <- tau
    
    ## sample from changepoint distribution given GP parameters
    prop <- as.numeric(rnorm(n = 1, mean = cp, sd = sqrt(cp_prop_var)))
    if(verbose == TRUE)
    {
      print(paste(i,"-th CP proposal: ", prop))
    }
    if(prop <= tol + interval[1] || prop >= -tol + interval[2])
    {
      par$cp[i + 1,] <- cp
    }
    else{
      temp_dat1 <- data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$y
      temp_dat2 <- data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$y
      
      prop_temp_dat1 <- data[data$x <= prop & data$x > interval[1], ]$y
      prop_temp_dat2 <- data[data$x < interval[2] & data$x > prop, ]$y
      
      log_accept_ratio <- lognormal_ou_pdf(x = prop_temp_dat1, mu = rep(0, times = length(prop_temp_dat1)), sigma = sigma[1], l = l[1]) +
        lognormal_ou_pdf(x = prop_temp_dat2, mu = rep(0, times = length(prop_temp_dat2)), sigma = sigma[2], l = l[2]) -
        (lognormal_ou_pdf(x = temp_dat1, mu = rep(0, times = length(temp_dat1)), sigma = sigma[1], l = l[1]) +
           lognormal_ou_pdf(x = temp_dat2, mu = rep(0, times = length(temp_dat2)), sigma = sigma[2], l = l[2]))
      if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
      {
        cp <- prop
      }
    }
    par$cp[i + 1,] <- cp
    #print(i)
  }
  ###########################################################
  ## End warmup
  ###########################################################
  ## tune metropolis proposal variances
  prop_var[[1]] <- 2.4^2 * var(cbind(par$sigma[warmup/2:warmup,1], par$l[warmup/2:warmup,1], par$tau[warmup/2:warmup,1])) / 3
  prop_var[[2]] <- 2.4^2 * var(cbind(par$sigma[warmup/2:warmup,2], par$l[warmup/2:warmup,2], par$tau[warmup/2:warmup,2])) / 3
  cp_prop_var <- 2.4^2 * var(par$cp[warmup/2:warmup,]) / 2
  
  ## reinitialize parameter list
  lp <- numeric()
  par <- list()
  par$sigma <- matrix(nrow = iter + 1, ncol = 2) ## the variance of the GP
  par$sigma[1,] <- sigma
  
  par$l <- matrix(nrow = iter + 1, ncol = 2) ## length scale of the GP
  par$l[1,] <- l
  
  par$tau <- matrix(nrow = iter + 1, ncol = 2) ## nugget of the data model
  par$tau[1,] <- tau
  
  par$cp <- matrix(nrow = iter + 1, ncol = 1) ## changepoint locations
  par$cp[1,] <- cp
  
  ## gibbs sampling iterations
  for(i in 1:(iter))
  {
    xrange <- matrix(nrow = 3, ncol = 2)
    xrange[1,] <- c(interval[1], cp[1])
    xrange[2,] <- c(cp[1], interval[2])
    
    ## given changepoints make proposal for MH steps for GP parameters
    for(j in 1:2)
    {
      prop <- as.numeric(rmvnorm(n = 1, mean = c(sigma[j], l[j], tau[j]), sigma = prop_var[[j]]))
      if(verbose == TRUE)
      {
        print(paste("iteration: ",i))
        print(paste(j,"-th GP parameter proposal: ", prop))
      }
      
      ## skip this chunk of data if the proposal produces negative values
      if(any(prop <= 0))
      {
        # par$sigma[i,j] <- sigma[j]
        # par$l[i,j] <- l[j]
        # par$tau[i,j] <- tau[j]
        next
      }
      else{
        temp_dat <- data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$y
        
        ## proposal doesn't appear because it should cancel
        log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = prop[1], l = prop[2]) +
          dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) +
          dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) -
          (lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma[j], l = l[j]) +
             dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) +
             dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE))
        
        if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
        {
          accept$gp_par[1,j] <- accept$gp_par[1,j] + 1/iter
          sigma[j] <- prop[1]
          l[j] <- prop[2]
          tau[j] <- prop[3]
        }
      }
    }
    ## update GP parameters
    par$sigma[i + 1,] <- sigma
    par$l[i + 1,] <- l
    par$tau[i + 1,] <- tau
    
    ## sample from changepoint distribution given GP parameters
    prop <- as.numeric(rnorm(n = 1, mean = cp, sd = sqrt(cp_prop_var)))
    if(verbose == TRUE)
    {
      print(paste(i,"-th CP proposal: ", prop))
    }
    if(prop[1] <= tol + interval[1] || prop[1] >= -tol + interval[2])
    {
      par$cp[i + 1,] <- cp
      temp_dat1 <- data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$y
      temp_dat2 <- data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$y
      lp[i] <- (lognormal_ou_pdf(x = temp_dat1, mu = rep(0, times = length(temp_dat1)), sigma = sigma[1], l = l[1]) +
                  lognormal_ou_pdf(x = temp_dat2, mu = rep(0, times = length(temp_dat2)), sigma = sigma[2], l = l[2]))
    }
    else{
      temp_dat1 <- data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$y
      temp_dat2 <- data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$y
      
      prop_temp_dat1 <- data[data$x <= prop & data$x > interval[1], ]$y
      prop_temp_dat2 <- data[data$x < interval[2] & data$x > prop, ]$y
      
      log_accept_ratio <- lognormal_ou_pdf(x = prop_temp_dat1, mu = rep(0, times = length(prop_temp_dat1)), sigma = sigma[1], l = l[1]) +
        lognormal_ou_pdf(x = prop_temp_dat2, mu = rep(0, times = length(prop_temp_dat2)), sigma = sigma[2], l = l[2]) -
        (lognormal_ou_pdf(x = temp_dat1, mu = rep(0, times = length(temp_dat1)), sigma = sigma[1], l = l[1]) +
           lognormal_ou_pdf(x = temp_dat2, mu = rep(0, times = length(temp_dat2)), sigma = sigma[2], l = l[2]))
      
      lp[i] <- (lognormal_ou_pdf(x = temp_dat1, mu = rep(0, times = length(temp_dat1)), sigma = sigma[1], l = l[1]) +
                  lognormal_ou_pdf(x = temp_dat2, mu = rep(0, times = length(temp_dat2)), sigma = sigma[2], l = l[2]))
      
      if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
      {
        cp <- prop
        accept$cp <- accept$cp + 1/iter
        lp[i] <- lognormal_ou_pdf(x = prop_temp_dat1, mu = rep(0, times = length(prop_temp_dat1)), sigma = sigma[1], l = l[1]) +
          lognormal_ou_pdf(x = prop_temp_dat2, mu = rep(0, times = length(prop_temp_dat2)), sigma = sigma[2], l = l[2])
      }
    }
    par$cp[i + 1,] <- cp
    #print(i)
  }
  
  return(list("parameters" = par, "accept" = accept, "lp" = lp))
}


