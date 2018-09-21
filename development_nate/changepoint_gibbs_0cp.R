cp0_gibbs <- function(data, iter, start.vals, prop_var, tol = 10, warmup = 5000, verbose = FALSE)
{
  ##data is a data frame with column x and column y
  
  ## initialize parameter list
  par <- list()
  par$sigma <- matrix(nrow = warmup + 1, ncol = 1) ## the variance of the GP
  par$sigma[1,] <- start.vals$sigma
  
  par$l <- matrix(nrow = warmup + 1, ncol = 1) ## length scale of the GP
  par$l[1,] <- start.vals$l
  
  par$tau <- matrix(nrow = warmup + 1, ncol = 1) ## nugget of the data model
  par$tau[1,] <- start.vals$tau
  
  ## range on the x-axis of data
  interval <- range(data$x)
  
  ## current values of parameters
  sigma <- start.vals$sigma
  l <- start.vals$l
  tau <- start.vals$tau
  
  ## initialize acceptance rates
  accept <- list()
  accept$gp_par <- matrix(data = c(0), nrow = 1, ncol = 1)
  
  ## gibbs warmup iterations
  for(i in 1:(warmup)){
    
  ## given changepoints make proposal for MH steps for GP parameters
  prop <- as.numeric(rmvnorm(n = 1, mean = c(sigma, l, tau), sigma = prop_var))
  if(verbose == TRUE)
  {
    print(paste("iteration: ",i))
    print(paste("GP parameter proposal: ", prop))
  }
  
  ## skip this chunk of data if the proposal produces negative values
  if(any(prop <= 0))
  {
    # par$sigma[i,j] <- sigma[j]
    # par$l[i,j] <- l[j]
    # par$tau[i,j] <- tau[j]
    
    ## update GP parameters
    par$sigma[i + 1,] <- sigma
    par$l[i + 1,] <- l
    par$tau[i + 1,] <- tau
  }
  else{
    temp_dat <- data$y
    
    ## proposal doesn't appear because it should cancel
    log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = prop[1], l = prop[2]) +
      dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) +
      dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) -
      (lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma, l = l) +
         dgamma(x = l, shape = 3, rate = 5, log = TRUE) +
         dnorm(x = sigma, mean = 0, sd = 1, log = TRUE))
    
    if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
    {
      sigma <- prop[1]
      l <- prop[2]
      tau <- prop[3]
    }
    ## update GP parameters
    par$sigma[i + 1,] <- sigma
    par$l[i + 1,] <- l
    par$tau[i + 1,] <- tau
  }
    #print(i)
  }
  ###########################################################
  ## End warmup
  ###########################################################
  ## tune metropolis proposal variances
  prop_var <- 2.4^2 * var(cbind(par$sigma[warmup/2:warmup,1], par$l[warmup/2:warmup,1], par$tau[warmup/2:warmup,1])) / 3
  
  ## reinitialize parameter list
  lp <- numeric()
  par <- list()
  par$sigma <- matrix(nrow = iter + 1, ncol = 1) ## the variance of the GP
  par$sigma[1,] <- sigma
  
  par$l <- matrix(nrow = iter + 1, ncol = 1) ## length scale of the GP
  par$l[1,] <- l
  
  par$tau <- matrix(nrow = iter + 1, ncol = 1) ## nugget of the data model
  par$tau[1,] <- tau
  
  
  ## gibbs sampling iterations
  for(i in 1:(iter))
  {
    ## given changepoints make proposal for MH steps for GP parameters
    prop <- as.numeric(rmvnorm(n = 1, mean = c(sigma, l, tau), sigma = prop_var))
    if(verbose == TRUE)
    {
      print(paste("iteration: ",i))
      print(paste("GP parameter proposal: ", prop))
    }
    
    ## skip this chunk of data if the proposal produces negative values
    if(any(prop <= 0))
    {
      ## update GP parameters
      par$sigma[i + 1,] <- sigma
      par$l[i + 1,] <- l
      par$tau[i + 1,] <- tau
      
      temp_dat <- data$y
      lp[i] <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma, l = l)
      
    }
    else{
      temp_dat <- data$y
      
      ## proposal doesn't appear because it should cancel
      log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = prop[1], l = prop[2]) +
        dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) +
        dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) -
        (lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma, l = l) +
           dgamma(x = l, shape = 3, rate = 5, log = TRUE) +
           dnorm(x = sigma, mean = 0, sd = 1, log = TRUE))
      
      lp[i] <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma, l = l)
      
      if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
      {
        accept$gp_par[1,] <- accept$gp_par[1,] + 1/iter
        sigma <- prop[1]
        l <- prop[2]
        tau <- prop[3]
        lp[i] <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = prop[1], l = prop[2])
      }
      ## update GP parameters
      par$sigma[i + 1,] <- sigma
      par$l[i + 1,] <- l
      par$tau[i + 1,] <- tau
    }
  }
  return(list("parameters" = par, "accept" = accept, "lp" = lp))
}
