mlgp <- function(y,x,tol = 1e-6)
{
  cov.fun <- function(x1, x2, par)
  {
    sigma <- exp(par[1])
    tau <- exp(par[2])
    return(sigma^2 * exp(-1/(2*tau^2) * (x1 - x2)^2))
  }
  make.cov.mat <- function(x, xpred = numeric(), cov.fun, par, eps = 0)
  {
    xfull <- c(x,xpred)
    temp <- matrix(nrow = length(xfull), ncol = length(xfull))
    for(i in 1:length(xfull))
    {
      for(j in 1:length(xfull))
      {
        temp[i,j] <- cov.fun(xfull[i], xfull[j], par = par) + 1*(i==j)*eps
      }
    }
    return(temp)
  }
  obj_fun <- function(par, ...)
  {
    ## make covariance matrix
    K <- make.cov.mat(x = x, xpred = numeric(), cov.fun = cov.fun, par = par, eps = 0)
    return(-dmvnorm(x = y, sigma = K, log = TRUE))
  }
  
  mle <- optim(par = c(0,1), fn = obj_fun, x = x, y = y)
  return(mle)
}

bullet_resid <- hamby44$ccdata_w_resid[[3]]
bullet_resid_ds <- bullet_resid[seq(from = 1, to = nrow(bullet_resid),by = 20),]
bullet_resid_ds <- bullet_resid_ds[complete.cases(bullet_resid_ds),]
test <- mlgp(y = bullet_resid_ds$rlo_resid, x = bullet_resid_ds$y)
test$par

myimpute <- function(y,x, sigma, l)
{
  cov.fun <- function(x1, x2, par)
  {
    sigma <- par$sigma
    tau <- par$tau
    return(sigma^2 * exp(-1/(2*tau^2) * (x1 - x2)^2))
  }
  
  make.cov.mat <- function(x, xpred, cov.fun, par, eps = 0)
  {
    xfull <- c(x,xpred)
    temp <- matrix(nrow = length(xfull), ncol = length(xfull))
    for(i in 1:length(xfull))
    {
      for(j in 1:length(xfull))
      {
        temp[i,j] <- cov.fun(xfull[i], xfull[j], par = par) + 1*(i==j)*eps
      }
    }
    return(temp)
  }
  normal.cond.mean <- function(y, x, xpred, mu = rep(0, times = length(c(x,xpred))), sigma)
  {
    sigma21 <- sigma[(length(x) + 1):length(mu),1:length(x)]
    sigma11 <- sigma[1:length(x), 1:length(x)]
    return(mu[(length(x) + 1):length(mu)] + sigma21 %*% qr.solve(sigma11) %*% (y - mu[1:length(x)]))
  }
  
  normal.cond.var <- function(y,x, xpred, sigma, mu)
  {
    sigma21 <- sigma[(length(x) + 1):length(mu),1:length(x)]
    sigma11 <- sigma[1:length(x), 1:length(x)]
    sigma22 <- sigma[(length(x) + 1):length(mu),(length(x) + 1):length(mu)]
    return(sigma22 - sigma21 %*% solve(sigma11) %*% t(sigma21))
  }
  
  ## split data into NA and non-NA data
  #y.NA <- y[is.na(y)]
  x.NA <- x[is.na(y)]
  y.ok <- y[!is.na(y)]
  x.ok <- x[!is.na(y)]
  x.ok <- x.ok[seq(from = 1, to = length(x.ok), by = 5)]
  y.ok <- y.ok[seq(from = 1, to = length(y.ok), by = 5)]
  
  ## make covariance matrix
  K <- make.cov.mat(x = x.ok, xpred = x.NA, cov.fun = cov.fun, par = list("sigma" = sigma, "tau" = l), eps = 1e-3)
  pred.y <- normal.cond.mean(y = y.ok, x = x.ok, xpred = x.NA, sigma = K)
  pred.df <- data.frame("x" = x, "y" = y)
  pred.df$y[is.na(y)] <- pred.y
  return(pred.df)
}

test.impute <- myimpute(y = bullet_resid$rlo_resid, x = bullet_resid$y, sigma = exp(test$par[1]), l = exp(test$par[2]))
plot(bullet_resid$y, bullet_resid$rlo_resid)
points(test.impute$x, test.impute$y, col = "red", type = "l")
