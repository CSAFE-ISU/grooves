## end to end data preprocessing and gibbs algorithm
detect_cp <- function(data, iter = 10000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol = 10, warmup = 5000, verbose = FALSE, prior_numcp = c(1/3, 1/3, 1/3), 
                      max_int_length = 10, max_cut_length = 200, int_length_interval = 1, cut_length_interval = 20)
{
  ## put extra functions in here just in case
  lognormal_ou_pdf <- function(x, mu, sigma, l)
  {
    n <- length(x)
    rho <- exp(-1/l)
    
    return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
           - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
                                                - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
  }
  avg_points <- function(y, x, windo)
  {
    avgs <- numeric()
    int_cntr <- numeric()
    ints <- seq(from = range(x)[1], to = range(x)[2], by = windo)
    for(i in 1:length(ints))
    {
      int_cntr[i] <- range(x)[1] + (i-1)*windo + windo/2
      avgs[i] <- mean(y[x > (range(x)[1] + (i-1)*windo) & x <= (range(x)[1] + (i)*windo)])
    }
    return(list("y" = avgs, "x" = int_cntr, "any.na" = any(is.na(avgs))))
  }
  
  preprocess_bullet_resid_simple <- function(data, max_int_length, max_cut_length, int_length_interval, cut_length_interval)
  {
    min_int_length <- min(data$x[2:nrow(data)] - data$x[1:(nrow(data) - 1)])
    min_cut_length <- 0
    
    int_length_vec <- seq(from = min_int_length, to = max_int_length, by = int_length_interval)
    cut_length_vec <- seq(from = min_cut_length, to = max_cut_length, by = cut_length_interval)
    na.array <- array(dim = c(length(int_length_vec), length(cut_length_vec), length(cut_length_vec)))
    obj.array <- array(dim = c(length(int_length_vec), length(cut_length_vec), length(cut_length_vec)))
    for(i in 1:length(int_length_vec))
    {
      for(j in 1:length(cut_length_vec))
      {
        for(k in 1:length(cut_length_vec))
        {
          obj.array[i,j,k] <- 1 - min(c(((range(data$x)[2] - range(data$x)[1] - cut_length_vec[j] - cut_length_vec[k])/int_length_vec[i]) / length(data$x),1))
          temp_range <- data$x >= min(data$x) + cut_length_vec[j] & data$x <= max(data$x) - cut_length_vec[k]
          na.array[i,j,k] <- (1e5)*avg_points(y = data[temp_range,]$y, x = data[temp_range,]$x, windo = int_length_vec[i])$any.na
        }
      }
    }
    ## Figure out which cut_length and interval_length produces the smallest error 
    ## while resulting in 0 NA's
    co <- which(x = (obj.array + na.array) == min(obj.array + na.array), arr.ind = TRUE)
    xrange <- c(range(data$x)[1] + cut_length_vec[co[1,2]], range(data$x)[2] - cut_length_vec[co[1,3]])
    d <- avg_points(y = data$y[data$x <= xrange[2] & data$x >= xrange[1]], 
                    x = data$x[data$x <= xrange[2] & data$x >= xrange[1]],
                    windo = int_length_vec[co[1,1]])
    
    return(list("obj.fun" = obj.array, "na" = na.array, "int_length" = int_length_vec[co[1,1]], "cut_length" = c(cut_length_vec[co[1,2]], cut_length_vec[co[1,3]]),
                "d" = d, "perc_observations" = nrow(d)/nrow(data), "perc_cut" = c(cut_length_vec[co[1,2]], cut_length_vec[co[1,3]])/range(data$x)))
  }
  ######### end extra functions
  
  ## remove NA values from the data 
  test_dat <- data[!is.na(data$rlo_resid),]
  d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
  
  ## preprocess data
  try(temp_d <- preprocess_bullet_resid_simple(data = d, max_int_length = max_int_length, max_cut_length = max_cut_length, int_length_interval = int_length_interval, cut_length_interval = cut_length_interval))

  ## run cp algorithm
  try(test_variable_cp_gibbs <- variable_cp_gibbs(data = data.frame("y" = temp_d$d$y, "x" = temp_d$d$x),
                                                          start.vals = start.vals, 
                                                          prop_var = prop_var, cp_prop_var = cp_prop_var, verbose = FALSE, tol = tol, iter = iter, warmup = warmup, prior_numcp = prior_numcp))
  return(list("changepoint_results" = test_variable_cp_gibbs, "cutoffs" = range(temp_d$d$x)))
}

# system.time(blob <- detect_cp(data = hamby44$ccdata_w_resid[[7]], warmup = 10000))
# plot(hamby44$ccdata_w_resid[[7]]$y, hamby44$ccdata_w_resid[[7]]$rlo_resid)
# blob$changepoint_results$posterior_numcp
# blob$changepoint_results$cp_mean
# blob$cutoffs
# abline(v = c(66,2200))
# abline(v = c(blob$changepoint_results$cp_mean$`2cp`), col = "red")
