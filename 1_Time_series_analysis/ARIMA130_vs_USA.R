library(rstan)
library(coda)
library(tseries)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)
library(readxl)
library(Metrics)
library(MLmetrics)

mydata <- read_excel("data.xlsx")
BIRTH_data <- data.frame(mydata)

#              Mean     SD Naive SE Time-series SE
# m0         0.1228 0.3502 0.007830       0.007827      good
# phi[1]    -0.8385 0.1377 0.003080       0.003122      good
# phi[2]    -0.5716 0.1828 0.004087       0.004058      good
# sigma2[1]  3.3809 1.0522 0.023527       0.025329      ok
# sigma2[2] 16.6686 5.1460 0.115067       0.109909      ok
set.seed(1441)
state <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", 
           "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia", "Hawaii", 
           "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky",
           "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
           "Mississippi", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey",
           "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma",
           "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
           "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia",
           "Wisconsin", "Wyoming")
N_try = 25
p <- c(0.1228, -0.8385, -0.5716, 3.3809, 16.6686)
mse_tot <- matrix(NA, nrow = length(state), ncol = 2)

for (s in 1:length(state)) {
  Yh_s <- BIRTH_data$State.Rate[BIRTH_data$State == state[s] & BIRTH_data$Age.Group..Years. == "15-17 years"]
  Yp_s <- BIRTH_data$State.Rate[BIRTH_data$State == state[s] & BIRTH_data$Age.Group..Years. == "18-19 years"]
  Y_real <- cbind(Yh_s, Yp_s)
  
  Y_1stdiff <- diff(Y_real,diff=1)
  Y_2nddiff <- diff(Y_real,diff=2)
  Y_3rddiff <- diff(Y_real,diff=3)
  
  Y_out_s <- Y_3rddiff[1,]
  res_h <- rep(NA, N_try)
  res_p <- rep(NA, N_try)
  
  for (i in 1:N_try) {
    if(i == 1){
      res_h[i] <- rnorm(1, p[1] + p[2] * Y_out_s[1], sqrt(p[4]))
      res_p[i] <- rnorm(1, p[1] + p[3] * Y_out_s[2], sqrt(p[5]))
    } else {
      res_h[i] <- rnorm(1, p[1] + p[2] * res_h[i-1], sqrt(p[4]))
      res_p[i] <- rnorm(1, p[1] + p[3] * res_p[i-1], sqrt(p[5]))
    }
  }
  
  Y_out_s <- rbind(Y_out_s, cbind(res_h,res_p))
  mse_h <- MSE(Y_out_s[,1], Y_3rddiff[,1])
  mse_p <- MSE(Y_out_s[,2], Y_3rddiff[,2])
  mse_tot[s,] <- c(mse_h, mse_p)
  
}

row.names(mse_tot) <- state
colnames(mse_tot) <- c("mse high", "mse post-high")
mse_tot




# we have to select some states from the mse table and then use these plots
ggplot(data.frame(x = 1:26), y = Y_out_s[, 1]) +
  geom_line(data = data.frame(Y = Y_out_s[, 1], x = 1:26), aes(x = x, y = Y), lty = 2, col = "red") +
  geom_line(data = data.frame(Y = Y_real_diff[, 1], x = 1:26), aes(x = x, y = Y), col = "blue") + theme_bw() +
  ggtitle("Forecasting on State ? for Highschool")

ggplot(data.frame(x = 1:26), y = Y_out_s[, 2]) +
  geom_line(data = data.frame(Y = Y_out_s[, 2], x = 1:26), aes(x = x, y = Y), lty = 2, col = "red") +
  geom_line(data = data.frame(Y = Y_real_diff[, 2], x = 1:26), aes(x = x, y = Y), col = "blue") + theme_bw() +
  ggtitle("Forecasting on State ? for Post-Highschool")


graphics.off()
rm(list=ls())
