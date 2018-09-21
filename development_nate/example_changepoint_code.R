## run the changepoint detection algorithm on the hamby set
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(locfit)
source(file = "development_nate/changepoint_gibbs_0cp.R")
source(file = "development_nate/changepoint_gibbs_1cp.R")
source(file = "development_nate/changepoint_gibbs_2cp.R")
source(file = "development_nate/variable_cp_gibbs.R")
source(file = "development_nate/end_to_end_cp_algorithm.R")


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

## add robust loess and robust lm residuals
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

## The following line of code runs the first version of the changepoint algorithm. 
## estimated changepoints are a vector with the name "grooves"
## This also takes about 2 minutes to run
test <- detect_cp(data = hamby44$ccdata_w_resid[[7]], warmup = 10000)

## groove estimates 
test$grooves

