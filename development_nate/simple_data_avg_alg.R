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

preprocess_bullet_resid_simple <- function(data, max_int_length = 10, max_cut_length = 200, int_length_interval = 1, cut_length_interval = 20)
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
         "d" = d, "perc_observations" = length(d$x)/nrow(data), "perc_cut" = c(cut_length_vec[co[1,2]], cut_length_vec[co[1,3]])/range(data$x)))
}


## hamby 1 
# bullet_resid <- hamby44$ccdata_w_resid[[1]]
# test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
# d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
# plot(d$x, d$y)
# system.time(test_preprocess <- preprocess_bullet_resid_simple(data = d))
# points(test_preprocess$d$x, test_preprocess$d$y, col = "red")
# test_preprocess$obj.fun
# test_preprocess$obj.fun + test_preprocess$na
# test_preprocess$int_length
# test_preprocess$cut_length
# 
# ## hamby 2
# bullet_resid <- hamby44$ccdata_w_resid[[2]]
# test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
# d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
# plot(d$x, d$y)
# system.time(test_preprocess <- preprocess_bullet_resid_simple(data = d))
# points(test_preprocess$d$x, test_preprocess$d$y, col = "red")
# test_preprocess$obj.fun
# test_preprocess$obj.fun + test_preprocess$na
# test_preprocess$int_length
# test_preprocess$cut_length
# 
# ## hamby 3
# bullet_resid <- hamby44$ccdata_w_resid[[3]]
# test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
# d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
# plot(d$x, d$y)
# system.time(test_preprocess <- preprocess_bullet_resid_simple(data = d))
# points(test_preprocess$d$x, test_preprocess$d$y, col = "red")
# test_preprocess$obj.fun
# test_preprocess$obj.fun + test_preprocess$na
# test_preprocess$int_length
# test_preprocess$cut_length
# 
# ## hamby 4
# bullet_resid <- hamby44$ccdata_w_resid[[4]]
# test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
# d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
# plot(d$x, d$y)
# system.time(test_preprocess <- preprocess_bullet_resid_simple(data = d))
# points(test_preprocess$d$x, test_preprocess$d$y, col = "red")
# test_preprocess$obj.fun
# test_preprocess$obj.fun + test_preprocess$na
# test_preprocess$int_length
# test_preprocess$cut_length
# 
# ## hamby 5
# bullet_resid <- hamby44$ccdata_w_resid[[5]]
# test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
# d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
# plot(d$x, d$y)
# system.time(test_preprocess <- preprocess_bullet_resid_simple(data = d))
# points(test_preprocess$d$x, test_preprocess$d$y, col = "red")
# test_preprocess$obj.fun
# test_preprocess$obj.fun + test_preprocess$na
# test_preprocess$int_length
# test_preprocess$cut_length
# 
## hamby 6
bullet_resid <- hamby44$ccdata_w_resid[[6]]
test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
plot(d$x, d$y)
system.time(test_preprocess <- preprocess_bullet_resid_simple(data = d))
points(test_preprocess$d$x, test_preprocess$d$y, col = "red")
test_preprocess$obj.fun
test_preprocess$obj.fun + test_preprocess$na
test_preprocess$int_length
test_preprocess$cut_length


