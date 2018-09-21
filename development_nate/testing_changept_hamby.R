## run the changepoint detection algorithm on the hamby set
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(locfit)

## create the robust residuals data sets
hamby44 <- readRDS("//opt/hamby44/hamby44.rda")

## filter NA values
hamby44 <- hamby44 %>% filter(!is.na(crosscuts))

## create a data frame with averaged crosscuts and subtract minimum value 
hamby44 <- hamby44 %>% mutate(ccdata_avg = purrr::map(ccdata, .f = function(dframe){
  dframe <- dframe %>% group_by(y) %>% summarise(value = mean(value,na.rm = TRUE))
  dframe <- as.data.frame(dframe)
  check_min <- min(dframe$value[!is.na(dframe$value)])
  dframe <- dframe %>% mutate(value_std = value - check_min)
  return(dframe)
}))

hamby44 <- hamby44 %>% mutate(ccdata_w_resid = purrr::map(ccdata_avg, .f = function(bullet){
  lm0 <- MASS::rlm(value_std~poly(y,2), data=bullet, maxit=100)
  bullet$pred <- predict(lm0, newdata=bullet)
  
  bullet$absresid <- with(bullet, abs(value_std-pred))
  bullet$resid <- with(bullet, value_std-pred)
  
  
  robust_loess_fit <- locfit.robust(value_std~y, data = bullet, alpha = 1, kern = "tcub")
  bullet$rlo_pred <- predict(robust_loess_fit, newdata = bullet)
  
  bullet$rlo_absresid <- with(bullet, abs(value_std-rlo_pred))
  bullet$rlo_resid <- with(bullet, value_std-rlo_pred)
  return(bullet)
}), rlm_r2 = purrr::map_dbl(ccdata_w_resid, .f = function(bullet){
  r2 <- (cor(bullet$value_std, bullet$pred, use = "complete.obs"))^2
  return(r2)
}), rlo_r2 = purrr::map_dbl(ccdata_w_resid, .f = function(bullet){
  r2 <- (cor(bullet$value_std, bullet$rlo_pred, use = "complete.obs"))^2
  return(r2)
})
)

## groove classes
hamby44$is_left_groove <- c(1, 1, 1, 1, 1, #5
                            1, 1, 1, 1, 0, #10
                            1, 1, 1, 1, 1, #15
                            1, 1, 1, 1, 1, #20
                            1, 1, 1, 1, 1, #25
                            1, 1, 1, 1, 1, #30
                            1, 0, 1, 1, 1, #35
                            1, 1, 1, 1, 1, #40
                            1, 1, 1, 1, 1, #45
                            1, 1, 1, 1, 1, #50
                            1, 1, 1, 1, 1, #55
                            1, 1, 0, 1, 1, #60
                            1, 1, 1, 1, 1, #65
                            1, 1, 1, 1, 1, #70
                            1, 1, 1, 1, 1, #75
                            1, 1, 1, 1, 1, #80
                            1, 1, 1, 1, 1, #85
                            1, 1, 1, 1, 1, #90
                            1, 1, 1, 1, 1, #95
                            1, 1, 1, 0, 1, #100
                            1, 1, 1, 1, 1, #105
                            1, 1, 1, 1, 1, #110
                            1, 1, 1, 1, 1, #115
                            1, 1, 1, 1, 1, #120
                            1, 1, 1, 1, 1, #125
                            1, 1, 1, 1, 1, #130
                            1, 1, 1, 1, 1, #135
                            1, 1, 1, 1, 1, #140
                            1, 1, 1, 1, 1, #145
                            1, 1, 1, 1, 1, #150
                            1, 1, 1, 1, 1, #155
                            1, 1, 1, 1, 1, #160
                            1, 1, 1, 1, 1, #165
                            1, 1, 1, 1, 1, #170
                            1, 1, 1, 1, 1, #175
                            1, 1, 1, 1, 1, #180
                            1, 1, 1, 1, 1, #185
                            1, 1, 1, 1, 1, #190
                            1, 1, 1, 1, 1, #195
                            1, 1, 1, 1, 1, #200
                            1, 1, 1, 1, 1, #205
                            1, 1, 1)

hamby44$is_right_groove <- c(1, 1, 0, 0, 0, #5 
                             1, 0, 0, 0, 1, #10
                             0, 0, 1, 1, 0, #15
                             0, 1, 1, 1, 0, #20
                             0, 0, 0, 1, 0, #25
                             1, 0, 0, 0, 0, #30
                             1, 1, 1, 0, 1, #35
                             0, 1, 1, 1, 0, #40
                             0, 0, 0, 1, 0, #45
                             1, 0, 1, 0, 0, #50
                             1, 0, 0, 1, 0, #55
                             0, 0, 1, 1, 0, #60
                             0, 1, 0, 0, 1, #65
                             1, 0, 1, 1, 1, #70
                             0, 0, 1, 0, 1, #75
                             0, 1, 0, 1, 0, #80
                             0, 1, 0, 1, 1, #85
                             1, 0, 1, 0, 1, #90
                             0, 0, 0, 1, 1, #95
                             1, 0, 1, 1, 1, #100
                             0, 0, 1, 1, 0, #105
                             0, 0, 1, 0, 0, #110
                             0, 1, 1, 1, 0, #115
                             0, 1, 0, 1, 0, #120
                             1, 1, 1, 0, 0, #125
                             0, 1, 1, 1, 1, #130
                             1, 1, 0, 1, 1, #135
                             1, 1, 0, 0, 0, #140
                             0, 0, 0, 0, 0, #145
                             0, 1, 1, 0, 0, #150
                             1, 1, 0, 1, 0, #155
                             0, 1, 0, 0, 0, #160
                             0, 0, 0, 0, 1, #165
                             1, 0, 1, 0, 0, #170
                             1, 1, 0, 0, 0, #175
                             0, 0, 1, 0, 1, #180
                             0, 1, 0, 1, 0, #185
                             1, 1, 1, 1, 0, #190
                             1, 1, 1, 0, 1, #195
                             1, 1, 1, 0, 1, #200
                             0, 1, 0, 1, 0, #205
                             1, 0, 1)

## run the changepoing algorithm and store the results in the purrr object(?)
system.time(gibberish <- hamby44 %>% mutate(changept_results = purrr::map(ccdata_w_resid, .f = detect_cp, warmup = 10000, tol = 40)))
plot(gibberish$ccdata_w_resid[[1]]$y, gibberish$ccdata_w_resid[[1]]$rlo_resid)
abline(v = gibberish$changept_results[[1]]$cutoffs)
abline(v = gibberish$changept_results[[1]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[1]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[1]]$changepoint_results$posterior_numcp
gibberish$changept_results[[1]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[2]]$y, gibberish$ccdata_w_resid[[2]]$rlo_resid)
abline(v = gibberish$changept_results[[2]]$cutoffs)
abline(v = gibberish$changept_results[[2]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[2]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[2]]$changepoint_results$posterior_numcp
gibberish$changept_results[[2]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[3]]$y, gibberish$ccdata_w_resid[[3]]$rlo_resid)
abline(v = gibberish$changept_results[[3]]$cutoffs)
abline(v = gibberish$changept_results[[3]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[3]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[3]]$changepoint_results$posterior_numcp
gibberish$changept_results[[3]]$changepoint_results$cp_mean

## bad result here 
## I suspect that the estimated variance in the middle segment is lowest
plot(gibberish$ccdata_w_resid[[4]]$y, gibberish$ccdata_w_resid[[4]]$rlo_resid)
abline(v = gibberish$changept_results[[4]]$cutoffs)
abline(v = gibberish$changept_results[[4]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[4]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[4]]$changepoint_results$posterior_numcp
gibberish$changept_results[[4]]$changepoint_results$cp_mean

## not great... why did the right cutoff move so far in?
plot(gibberish$ccdata_w_resid[[5]]$y, gibberish$ccdata_w_resid[[5]]$rlo_resid)
abline(v = gibberish$changept_results[[5]]$cutoffs)
abline(v = gibberish$changept_results[[5]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[5]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[5]]$changepoint_results$posterior_numcp
gibberish$changept_results[[5]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[6]]$y, gibberish$ccdata_w_resid[[6]]$rlo_resid)
abline(v = gibberish$changept_results[[6]]$cutoffs)
abline(v = gibberish$changept_results[[6]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[6]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[6]]$changepoint_results$posterior_numcp
gibberish$changept_results[[6]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[7]]$y, gibberish$ccdata_w_resid[[7]]$rlo_resid)
abline(v = gibberish$changept_results[[7]]$cutoffs)
abline(v = gibberish$changept_results[[7]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[7]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[7]]$changepoint_results$posterior_numcp
gibberish$changept_results[[7]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[8]]$y, gibberish$ccdata_w_resid[[8]]$rlo_resid)
abline(v = gibberish$changept_results[[8]]$cutoffs)
abline(v = gibberish$changept_results[[8]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[8]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[8]]$changepoint_results$posterior_numcp
gibberish$changept_results[[8]]$changepoint_results$cp_mean

## bad result here
## the two changepoints should not have been allowed to be so close to each other
plot(gibberish$ccdata_w_resid[[9]]$y, gibberish$ccdata_w_resid[[9]]$rlo_resid)
abline(v = gibberish$changept_results[[9]]$cutoffs)
abline(v = gibberish$changept_results[[9]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[9]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[9]]$changepoint_results$posterior_numcp
gibberish$changept_results[[9]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[10]]$y, gibberish$ccdata_w_resid[[10]]$rlo_resid)
abline(v = gibberish$changept_results[[10]]$cutoffs)
abline(v = gibberish$changept_results[[10]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[10]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[10]]$changepoint_results$posterior_numcp
gibberish$changept_results[[10]]$changepoint_results$cp_mean

## why does it think that there are two changepoints here?
plot(gibberish$ccdata_w_resid[[11]]$y, gibberish$ccdata_w_resid[[11]]$rlo_resid)
abline(v = gibberish$changept_results[[11]]$cutoffs)
abline(v = gibberish$changept_results[[11]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[11]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[11]]$changepoint_results$posterior_numcp
gibberish$changept_results[[11]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[12]]$y, gibberish$ccdata_w_resid[[12]]$rlo_resid)
abline(v = gibberish$changept_results[[12]]$cutoffs)
abline(v = gibberish$changept_results[[12]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[12]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[12]]$changepoint_results$posterior_numcp
gibberish$changept_results[[12]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[13]]$y, gibberish$ccdata_w_resid[[13]]$rlo_resid)
abline(v = gibberish$changept_results[[13]]$cutoffs)
abline(v = gibberish$changept_results[[13]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[13]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[13]]$changepoint_results$posterior_numcp
gibberish$changept_results[[13]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[14]]$y, gibberish$ccdata_w_resid[[14]]$rlo_resid)
abline(v = gibberish$changept_results[[14]]$cutoffs)
abline(v = gibberish$changept_results[[14]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[14]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[14]]$changepoint_results$posterior_numcp
gibberish$changept_results[[14]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[15]]$y, gibberish$ccdata_w_resid[[15]]$rlo_resid)
abline(v = gibberish$changept_results[[15]]$cutoffs)
abline(v = gibberish$changept_results[[15]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[15]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[15]]$changepoint_results$posterior_numcp
gibberish$changept_results[[15]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[16]]$y, gibberish$ccdata_w_resid[[16]]$rlo_resid)
abline(v = gibberish$changept_results[[16]]$cutoffs)
abline(v = gibberish$changept_results[[16]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[16]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[16]]$changepoint_results$posterior_numcp
gibberish$changept_results[[16]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[17]]$y, gibberish$ccdata_w_resid[[17]]$rlo_resid)
abline(v = gibberish$changept_results[[17]]$cutoffs)
abline(v = gibberish$changept_results[[17]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[17]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[17]]$changepoint_results$posterior_numcp
gibberish$changept_results[[17]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[18]]$y, gibberish$ccdata_w_resid[[18]]$rlo_resid)
abline(v = gibberish$changept_results[[18]]$cutoffs)
abline(v = gibberish$changept_results[[18]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[18]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[18]]$changepoint_results$posterior_numcp
gibberish$changept_results[[18]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[19]]$y, gibberish$ccdata_w_resid[[19]]$rlo_resid)
abline(v = gibberish$changept_results[[19]]$cutoffs)
abline(v = gibberish$changept_results[[19]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[19]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[19]]$changepoint_results$posterior_numcp
gibberish$changept_results[[19]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[20]]$y, gibberish$ccdata_w_resid[[20]]$rlo_resid)
abline(v = gibberish$changept_results[[20]]$cutoffs)
abline(v = gibberish$changept_results[[20]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[20]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[20]]$changepoint_results$posterior_numcp
gibberish$changept_results[[20]]$changepoint_results$cp_mean

## bad
plot(gibberish$ccdata_w_resid[[21]]$y, gibberish$ccdata_w_resid[[21]]$rlo_resid)
abline(v = gibberish$changept_results[[21]]$cutoffs)
abline(v = gibberish$changept_results[[21]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[21]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[21]]$changepoint_results$posterior_numcp
gibberish$changept_results[[21]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[22]]$y, gibberish$ccdata_w_resid[[22]]$rlo_resid)
abline(v = gibberish$changept_results[[22]]$cutoffs)
abline(v = gibberish$changept_results[[22]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[22]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[22]]$changepoint_results$posterior_numcp
gibberish$changept_results[[22]]$changepoint_results$cp_mean

## really bad! Interesting case here.
plot(gibberish$ccdata_w_resid[[23]]$y, gibberish$ccdata_w_resid[[23]]$rlo_resid)
abline(v = gibberish$changept_results[[23]]$cutoffs)
abline(v = gibberish$changept_results[[23]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[23]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[23]]$changepoint_results$posterior_numcp
gibberish$changept_results[[23]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[24]]$y, gibberish$ccdata_w_resid[[24]]$rlo_resid)
abline(v = gibberish$changept_results[[24]]$cutoffs)
abline(v = gibberish$changept_results[[24]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[24]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[24]]$changepoint_results$posterior_numcp
gibberish$changept_results[[24]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[25]]$y, gibberish$ccdata_w_resid[[25]]$rlo_resid)
abline(v = gibberish$changept_results[[25]]$cutoffs)
abline(v = gibberish$changept_results[[25]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[25]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[25]]$changepoint_results$posterior_numcp
gibberish$changept_results[[25]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[26]]$y, gibberish$ccdata_w_resid[[26]]$rlo_resid)
abline(v = gibberish$changept_results[[26]]$cutoffs)
abline(v = gibberish$changept_results[[26]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[26]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[26]]$changepoint_results$posterior_numcp
gibberish$changept_results[[26]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[27]]$y, gibberish$ccdata_w_resid[[27]]$rlo_resid)
abline(v = gibberish$changept_results[[27]]$cutoffs)
abline(v = gibberish$changept_results[[27]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[27]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[27]]$changepoint_results$posterior_numcp
gibberish$changept_results[[27]]$changepoint_results$cp_mean

## nice 
plot(gibberish$ccdata_w_resid[[28]]$y, gibberish$ccdata_w_resid[[28]]$rlo_resid)
abline(v = gibberish$changept_results[[28]]$cutoffs)
abline(v = gibberish$changept_results[[28]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[28]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[28]]$changepoint_results$posterior_numcp
gibberish$changept_results[[28]]$changepoint_results$cp_mean

## bad case
plot(gibberish$ccdata_w_resid[[29]]$y, gibberish$ccdata_w_resid[[29]]$rlo_resid)
abline(v = gibberish$changept_results[[29]]$cutoffs)
abline(v = gibberish$changept_results[[29]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[29]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[29]]$changepoint_results$posterior_numcp
gibberish$changept_results[[29]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[30]]$y, gibberish$ccdata_w_resid[[30]]$rlo_resid)
abline(v = gibberish$changept_results[[30]]$cutoffs)
abline(v = gibberish$changept_results[[30]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[30]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[30]]$changepoint_results$posterior_numcp
gibberish$changept_results[[30]]$changepoint_results$cp_mean

## bad case
plot(gibberish$ccdata_w_resid[[31]]$y, gibberish$ccdata_w_resid[[31]]$rlo_resid)
abline(v = gibberish$changept_results[[31]]$cutoffs)
abline(v = gibberish$changept_results[[31]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[31]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[31]]$changepoint_results$posterior_numcp
gibberish$changept_results[[31]]$changepoint_results$cp_mean

## bad
plot(gibberish$ccdata_w_resid[[32]]$y, gibberish$ccdata_w_resid[[32]]$rlo_resid)
abline(v = gibberish$changept_results[[32]]$cutoffs)
abline(v = gibberish$changept_results[[32]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[32]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[32]]$changepoint_results$posterior_numcp
gibberish$changept_results[[32]]$changepoint_results$cp_mean


plot(gibberish$ccdata_w_resid[[33]]$y, gibberish$ccdata_w_resid[[33]]$rlo_resid)
abline(v = gibberish$changept_results[[33]]$cutoffs)
abline(v = gibberish$changept_results[[33]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[33]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[33]]$changepoint_results$posterior_numcp
gibberish$changept_results[[33]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[34]]$y, gibberish$ccdata_w_resid[[34]]$rlo_resid)
abline(v = gibberish$changept_results[[34]]$cutoffs)
abline(v = gibberish$changept_results[[34]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[34]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[34]]$changepoint_results$posterior_numcp
gibberish$changept_results[[34]]$changepoint_results$cp_mean

## bad
plot(gibberish$ccdata_w_resid[[35]]$y, gibberish$ccdata_w_resid[[35]]$rlo_resid)
abline(v = gibberish$changept_results[[35]]$cutoffs)
abline(v = gibberish$changept_results[[35]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[35]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[35]]$changepoint_results$posterior_numcp
gibberish$changept_results[[35]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[36]]$y, gibberish$ccdata_w_resid[[36]]$rlo_resid)
abline(v = gibberish$changept_results[[36]]$cutoffs)
abline(v = gibberish$changept_results[[36]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[36]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[36]]$changepoint_results$posterior_numcp
gibberish$changept_results[[36]]$changepoint_results$cp_mean

## really bad
plot(gibberish$ccdata_w_resid[[37]]$y, gibberish$ccdata_w_resid[[37]]$rlo_resid)
abline(v = gibberish$changept_results[[37]]$cutoffs)
abline(v = gibberish$changept_results[[37]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[37]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[37]]$changepoint_results$posterior_numcp
gibberish$changept_results[[37]]$changepoint_results$cp_mean

## bad
plot(gibberish$ccdata_w_resid[[38]]$y, gibberish$ccdata_w_resid[[38]]$rlo_resid)
abline(v = gibberish$changept_results[[38]]$cutoffs)
abline(v = gibberish$changept_results[[38]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[38]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[38]]$changepoint_results$posterior_numcp
gibberish$changept_results[[38]]$changepoint_results$cp_mean

## bad
plot(gibberish$ccdata_w_resid[[39]]$y, gibberish$ccdata_w_resid[[39]]$rlo_resid)
abline(v = gibberish$changept_results[[39]]$cutoffs)
abline(v = gibberish$changept_results[[39]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[39]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[39]]$changepoint_results$posterior_numcp
gibberish$changept_results[[39]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[40]]$y, gibberish$ccdata_w_resid[[40]]$rlo_resid)
abline(v = gibberish$changept_results[[40]]$cutoffs)
abline(v = gibberish$changept_results[[40]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[40]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[40]]$changepoint_results$posterior_numcp
gibberish$changept_results[[40]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[41]]$y, gibberish$ccdata_w_resid[[41]]$rlo_resid)
abline(v = gibberish$changept_results[[41]]$cutoffs)
abline(v = gibberish$changept_results[[41]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[41]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[41]]$changepoint_results$posterior_numcp
gibberish$changept_results[[41]]$changepoint_results$cp_mean

plot(gibberish$ccdata_w_resid[[41]]$y, gibberish$ccdata_w_resid[[41]]$rlo_resid)
abline(v = gibberish$changept_results[[41]]$cutoffs)
abline(v = gibberish$changept_results[[41]]$changepoint_results$cp_mean[[2]], col = "red")
abline(v = gibberish$changept_results[[41]]$changepoint_results$cp_mean[[1]], col = "blue")
gibberish$changept_results[[41]]$changepoint_results$posterior_numcp
gibberish$changept_results[[41]]$changepoint_results$cp_mean

