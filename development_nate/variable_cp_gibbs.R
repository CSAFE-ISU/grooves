lognormal_ou_pdf <- function(x, mu, sigma, l)
{
  n <- length(x)
  rho <- exp(-1/l)
  
  return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
         - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
                                              - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
}
lognormal_ou_pdf_par <- function(par, x)
{
  sigma <- par[1]
  l <- par[2]
  mu <- rep(0, times = length(x))
  n <- length(x)
  rho <- exp(-1/l)
  
  return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
         - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
                                              - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
}

## function to sample data weighted according to distance from center
sample_weights <- function(x, sq = FALSE, offset = 1)
{
  center <- (min(x) + max(x))/2
  ifelse(sq == FALSE, dist <- abs(x - center), dist <- (abs(x - center)^2))
  weights <- (dist + offset)/sum(dist + offset)
  return(weights)
}

## function to get the conditional posterior given 0,1,2 changepoints
variable_cp_gibbs <- function(data, iter = 10000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol = 10, warmup = 5000, verbose = FALSE, prior_numcp = c(1/3, 1/3, 1/3))
{
  ## If some function arguments (starting values/proposal variances are unspecified)
  ## choose generic arguments
  if(is.na(start.vals))
  {
    cp_start <- sort(sample(x = seq(from = range(data$x)[1], to = range(data$x)[2], by = 1), size = 2, replace = FALSE, prob = sample_weights(x = seq(from = range(data$x)[1], to = range(data$x)[2], by = 1))))
    start.vals <- list("cp2" = list("sigma" = c(1,1,1), "l" = c(10,10,10), "tau" = c(1,1,1), "cp" = cp_start),
                       "cp1" = list("sigma" = c(1,1), "l" = c(10,10), "tau" = c(1,1), "cp" = c(1000)),
                       "cp0" = list("sigma" = c(1), "l" = c(10), "tau" = c(1)))
  }
  if(is.na(prop_var))
  {
    prop_var <- list("cp2" = list(diag(c(1,1,1)), diag(c(1,1,1)), diag(c(1,1,1))),
                     "cp1" = list(diag(c(1,1,1)), diag(c(1,1,1))),
                     "cp0" = diag(c(1,1,1)))
  }
  if(is.na(cp_prop_var))
  {
    cp_prop_var <- list("cp2" = diag(c(400^2, 400^2)),
                        "cp1" = 400^2)
  }
  
  ## change point parameter list
  cp_list <- list()
  
  ## two changepoint model
  cp2_dsn <- cp2_gibbs(data = data, iter = iter, start.vals = start.vals$cp2, prop_var = prop_var$cp2, cp_prop_var = cp_prop_var$cp2, tol = tol, warmup = warmup, verbose = verbose)
  cp_list$cp2 <- cp2_dsn$parameters$cp
  mcp2 <- mean(cp2_dsn$lp)
  
  ## one changepoint model
  cp1_dsn <- cp1_gibbs(data = data, iter = iter, start.vals = start.vals$cp1, prop_var = prop_var$cp1, cp_prop_var = cp_prop_var$cp1, tol = tol, warmup = warmup, verbose = verbose)
  cp_list$cp1 <- cp1_dsn$parameters$cp
  mcp1 <- mean(cp1_dsn$lp)
  
  ## zero changepoint model
  cp0_dsn <- cp0_gibbs(data = data, iter = iter, start.vals = start.vals$cp0, prop_var = prop_var$cp0, tol = tol, warmup = warmup, verbose = verbose)
  mcp0 <- mean(cp0_dsn$lp)
  
  ## posterior cp probabilities
  ratio02 <- exp(log(prior_numcp[1]) + mcp0 - log(prior_numcp[3]) - mcp2)
  ratio12 <- exp(log(prior_numcp[2]) + mcp1 - log(prior_numcp[3]) - mcp2)
  
  p2 <- 1/(ratio02 + ratio12 + 1)
  p1 <- ratio12 * p2
  p0 <- ratio02 * p2
  
  post_numcp <- c(p0,p1,p2)
  
  return(list("posterior_numcp" = post_numcp, "posterior_cp" = cp_list,
              "cp_mean" = list("2cp" = apply(X = cp_list$cp2, MARGIN = 2, FUN = mean),
                               "1cp" = mean(cp_list$cp1))))
    
}

# ## test the full procedure
# bullet_resid <- hamby44$ccdata_w_resid[[6]]
# test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
# d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
# plot(d$x, d$y)
# 
# ## preprocess data
# temp_d <- preprocess_bullet_resid_simple(data = d)
# points(temp_d$d$x, temp_d$d$y, col = "red")
# 
# ## run cp algorithm
# system.time(test_variable_cp_gibbs <- variable_cp_gibbs(data = data.frame("y" = temp_d$d$y, "x" = temp_d$d$x),
#                                             start.vals = start.vals, 
#                                             prop_var = prop_var, cp_prop_var = cp_prop_var, verbose = FALSE))
# test_variable_cp_gibbs$cp_mean
# test_variable_cp_gibbs$posterior_numcp
