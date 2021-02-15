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
set.seed(14)

diff_to_timeseries <- function(y_last, y_1stdiff_pred) {
  a<- y_last
  y_dix_pred <-rep(0,length(y_1stdiff_pred))
  for (i in 1:length(y_1stdiff_pred)){
    y_dix_pred[i] <- a + y_1stdiff_pred[i]
    a <- y_dix_pred[i] 
  }
  return(y_dix_pred)
}

# Dataset construction
mydata <- read_excel("data.xlsx")
BIRTH_data <- data.frame(mydata)
Yh <- BIRTH_data$State.Rate[BIRTH_data$State == "Missouri" & BIRTH_data$Age.Group..Years. == "15-17 years"]
Yp <- BIRTH_data$State.Rate[BIRTH_data$State == "Missouri" & BIRTH_data$Age.Group..Years. == "18-19 years"]
Y <- cbind(Yh,Yp)
N <- dim(Y)[1]
A <- dim(Y)[2]

# Times series plot
obs <- seq(1990,2018,1)
y <- data.frame(Y)
ggplot(y) + geom_line(aes(x=obs, y=Y[,1], color="high schooler")) + geom_line(aes(x=obs, y=Y[,2], col="post high schooler")) + scale_color_discrete(name="Legend") + labs(title="Time Series")

# ADF test:
Y_1stdiff <- diff(Y,diff=1)
adf.test(Y_1stdiff[,1]) # 0.284
adf.test(Y_1stdiff[,2]) # 0.2237
Y_2nddiff <- diff(Y,diff=2)
adf.test(Y_2nddiff[,1]) # 0.201
adf.test(Y_2nddiff[,2]) # 0.09236
Y_3rddiff <- diff(Y,diff=3)
adf.test(Y_3rddiff[,1]) # 0.01
adf.test(Y_3rddiff[,2]) # 0.01116

# Our choice
Y_selected <- Y_3rddiff

# Training and Test set
N_test = 3
N_train = dim(Y_selected)[1] - N_test 
Y_train <- Y_selected[1:N_train,]# we work with 3rd order differences, so we have 23 obs
Y_test <- Y_selected[-(1:N_train),]

# plot of training set
obs <- obs[1:N_train]
y_train <- data.frame(Y_train)
ggplot(y_train) + geom_line(aes(x=obs, y=Y_train[,1], color="high schooler")) + geom_line(aes(x=obs, y=Y_train[,2], col="post high schooler")) + scale_color_discrete(name="Legend") + labs(title="Time Series of differences")


# ACF & PACF 
par(mfrow=c(1,2))
acf(Y_train[,1],lag =20)
pacf(Y_train[,1],lag =20)
# it suggests ARMA(2,3)

par(mfrow=c(1,2))
acf(Y_train[,2],lag =20)
pacf(Y_train[,2],lag =20)
# it suggests AR(1)

########################### FIRST MODEL: ARIMA(1,3,0) ############################
#------------------------------------------
# y_tj = alpha + phi_j * y_t-1,j + eps_j        AR(1) for 3rd difference data

# creating the list and evaluating the model
data_StSp <-list(N = N_train,
                 A = A,
                 Y = Y_train,
                 sigma2phi1 = 0.5, 
                 sigma2phi2 = 0.5, 
                 a_sigma2 = 2,
                 b_sigma2 = 10,
                 sigma2m0 = 10)

inits <- function() 
{
  list(phi = c(0,0),
       m0 = 0, 
       sigma2 = c(var(Yh),var(Yp)))
}

# StSp_model <- rstan::stan(file = "StSp_model_2ages.stan",
#                    data = data_StSp,
#                    chains = 2,
#                    iter = 30000,
#                    warmup = 10000,
#                    thin= 20,
#                    seed = 42,
#                    init = inits,
#                    algorithm = 'NUTS')
# 
# save(StSp_model, file="AR1_mixed_ages_3rd_diff.dat")
load("AR1_mixed_ages_3rd_diff.dat")

#x11()
rstan::traceplot(StSp_model, pars = c("m0", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"), inc_warmup = TRUE)


# Diagnostic

# coda
coda_chain <- rstan::As.mcmc.list(StSp_model, pars = c("m0", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"))
summary(coda_chain)
#              Mean     SD Naive SE Time-series SE
# m0         0.1228 0.3502 0.007830       0.007827      good
# phi[1]    -0.8385 0.1377 0.003080       0.003122      good
# phi[2]    -0.5716 0.1828 0.004087       0.004058      good
# sigma2[1]  3.3809 1.0522 0.023527       0.025329      ok
# sigma2[2] 16.6686 5.1460 0.115067       0.109909      ok

#              2.5%     25%     50%     75%   97.5%
# m0        -0.5552 -0.1223  0.1131  0.3685  0.8329
# phi[1]    -1.1213 -0.9240 -0.8393 -0.7495 -0.5694     NOT true that |phi| < 1 in 95% CI
# phi[2]    -0.9215 -0.7010 -0.5739 -0.4484 -0.2083     true that |phi| < 1 in 95% CI
# sigma2[1]  1.8719  2.6344  3.1823  3.9294  5.8973
# sigma2[2]  9.4872 13.0719 15.7557 19.2179 29.2515

# gelman and rubin's convergence diagnostic
coda::gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
coda::geweke.diag(coda_chain, frac1=0.1, frac2=0.5)
# ok, for parameters ~ Normal we have 1

# Autocorrelation and plot
#x11()
rstan::stan_ac(StSp_model, pars = c("m0", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"))
# ok, seems we can assume independent samples


# Posterior distributions
plot_post <- StSp_model %>% 
  rstan::extract(c("m0", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

#x11()
plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# Predictions
params <- data.frame(rstan::extract(StSp_model, c("m0", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"), perm = T))
resultStSp_1 <- matrix(NA, nrow = nrow(params), ncol = N_test)
resultStSp_2 <- matrix(NA, nrow = nrow(params), ncol = N_test)
for (i in 1:N_test) {
  if(i == 1){
    resultStSp_1[,i] <- apply(params, 1, function(x) rnorm(1, x[1] + x[2] * Y_train[N_train,1], sqrt(x[4])))
    resultStSp_2[,i] <- apply(params, 1, function(x) rnorm(1, x[1] + x[3] * Y_train[N_train,2], sqrt(x[5])))
  } else {
    resultStSp_1[,i] <- apply(cbind(params, resultStSp_1[,i-1]), 1, function(x) rnorm(1, x[1] + x[2] * x[6], sqrt(x[4])))
    resultStSp_2[,i] <- apply(cbind(params, resultStSp_2[,i-1]), 1, function(x) rnorm(1, x[1] + x[3] * x[6], sqrt(x[5])))
  }
}

ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,1],colMeans(resultStSp_1)),
                  ylow = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,1], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,1],Y_test[,1]) )) + 
  ggtitle("AR1 model - Third order differences data Missouri, mixed ages, high schooler")

ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,2],colMeans(resultStSp_2)),
                  ylow = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,2], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,2],Y_test[,2]) )) + 
  ggtitle("AR1 model - Third order differences data Missouri, mixed ages, post high schooler")


Y_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],colMeans(resultStSp_1))
Y_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_2nddiff_pred_1)
Y_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_1stdiff_pred_1)

Y_low_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_low_2nddiff_pred_1)
Y_low_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_low_1stdiff_pred_1)

Y_up_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_up_2nddiff_pred_1)
Y_up_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_up_1stdiff_pred_1)

#x11()
ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,1],Y_pred_1), 
                  ylow = c(Y[26,1],Y_low_pred_1), 
                  yup = c(Y[26,1],Y_up_pred_1)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Yh[1:26], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,1] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, high schooler")



Y_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],colMeans(resultStSp_2))
Y_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_2nddiff_pred_2)
Y_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_1stdiff_pred_2)

Y_low_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_low_2nddiff_pred_2)
Y_low_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_low_1stdiff_pred_2)

Y_up_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_up_2nddiff_pred_2)
Y_up_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_up_1stdiff_pred_2)

#x11()
ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,2],Y_pred_2), 
                  ylow = c(Y[26,2],Y_low_pred_2), 
                  yup = c(Y[26,2],Y_up_pred_2)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y[1:26,2], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,2] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, post high schooler")


# Goodness of fit
llik   <- data.frame(rstan::extract(StSp_model, "log_lik")[[1]])

p_WAIC_1 <- sum(apply(llik[,1:23], 2, var))
lppd_1   <- sum(apply(llik[,1:23], 2, function(x) log(mean(exp(x)))))
WAIC_1   <- - 2 * lppd_1 + 2 * p_WAIC_1   ##91.25328

p_WAIC_2 <- sum(apply(llik[,24:46], 2, var))
lppd_2   <- sum(apply(llik[,24:46], 2, function(x) log(mean(exp(x)))))
WAIC_2   <- - 2 * lppd_2 + 2 * p_WAIC_2   ##134.1255

# MSE for evaluate forcasting performance
MSE_1 <- MSE(Y_pred_1, Y[27:29,1])  #1.823686
MSE_2 <- MSE(Y_pred_2, Y[27:29,2])  #1.528984
tab <- rbind(c(WAIC_1, MSE_1, WAIC_2, MSE_2))
#------------------------------------ END ------------------------------------#



##################################################################################


########################### SECOND MODEL: ARIMA(1,3,0) ############################
#------------------------------------------
# y_tj = alpha_j + phi * y_t-1,j + eps_j        AR(1) for 3rd difference data

data_StSp2 <-list(N = N_train,
                  A = A,
                  Y = Y_train,
                  sigma2m0_1 = 10,
                  sigma2m0_2 = 10,
                  sigma2phi = 10,
                  a_sigma2 = 2,
                  b_sigma2 = 10)

inits <- function() 
{
  list(phi = 0,
       m0 = c(0,0), 
       sigma2 = c(var(Yh),var(Yp)))
}

# StSp_model2 <- rstan::stan(file = "StSp_model_2ages_second_version.stan",
#                    data = data_StS2,
#                    chains = 2,
#                    iter = 30000,
#                    warmup = 10000,
#                    thin= 20,
#                    seed = 42,
#                    init = inits,
#                    algorithm = 'NUTS')
# 
# save(StSp_model, file="AR1_mixed_ages_3rd_diff_2nd_version.dat")
load("AR1_mixed_ages_3rd_diff_2nd_version.dat")

rstan::traceplot(StSp_model2, pars = c("m0[1]", "m0[2]", "phi", "sigma2[1]", "sigma2[2]"), inc_warmup = TRUE)

# Diagnostic 

# coda
coda_chain <- rstan::As.mcmc.list(StSp_model2, pars = c("m0[1]", "m0[2]", "phi", "sigma2[1]", "sigma2[2]"))
summary(coda_chain)
#               Mean     SD Naive SE Time-series SE
# m0[1]      0.04698 0.3879 0.008674       0.008563       good
# m0[2]      0.57434 0.8960 0.020035       0.020040       good
# phi       -0.74758 0.1116 0.002494       0.002431       good
# sigma2[1]  3.40321 1.0764 0.024069       0.026139       ok
# sigma2[2] 17.21190 5.6075 0.125388       0.119130       ok

#              2.5%      25%      50%     75%   97.5%
# m0[1]     -0.7046 -0.19877  0.04189  0.3085  0.7945
# m0[2]     -1.1406 -0.02991  0.54293  1.1712  2.3445
# phi       -0.9678 -0.81709 -0.74715 -0.6733 -0.5322   true that |phi| < 1 in a 95% CI
# sigma2[1]  1.9225  2.64648  3.17874  3.9524  6.0576
# sigma2[2]  9.6239 13.38504 16.12691 19.8661 30.2097

# Gelman and rubin's convergence diagnostic
coda::gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
coda::geweke.diag(coda_chain, frac1=0.1, frac2=0.5)
# ok, parameters ~ Normal have 1 as an estimate

# Autocorrelation and plot
rstan::stan_ac(StSp_model2, pars = c("m0[1]", "m0[2]", "phi", "sigma2[1]", "sigma2[2]"))
# ok, we have independent samples

# Posterior distributions
plot_post <- StSp_model2 %>% 
  rstan::extract(c("m0[1]", "m0[2]", "phi", "sigma2[1]", "sigma2[2]")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")


# Predictions
params <- data.frame(rstan::extract(StSp_model2, c("m0[1]", "m0[2]", "phi", "sigma2[1]", "sigma2[2]"), perm = T))
resultStSp_1 <- matrix(NA, nrow = nrow(params), ncol = N_test)
resultStSp_2 <- matrix(NA, nrow = nrow(params), ncol = N_test)
for (i in 1:N_test) {
  if(i == 1){
    resultStSp_1[,i] <- apply(params, 1, function(x) rnorm(1, x[1] + x[3] * Y_train[N_train,1], sqrt(x[4])))
    resultStSp_2[,i] <- apply(params, 1, function(x) rnorm(1, x[2] + x[3] * Y_train[N_train,2], sqrt(x[5])))
  } else {
    resultStSp_1[,i] <- apply(cbind(params, resultStSp_1[,i-1]), 1, function(x) rnorm(1, x[1] + x[3] * x[6], sqrt(x[4])))
    resultStSp_2[,i] <- apply(cbind(params, resultStSp_2[,i-1]), 1, function(x) rnorm(1, x[2] + x[3] * x[6], sqrt(x[5])))
  }
}

ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,1],colMeans(resultStSp_1)),
                  ylow = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,1], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,1],Y_test[,1]) )) + 
  ggtitle("AR1 model - 3rd order differences data Missouri, mixed ages, high schooler")

ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,2],colMeans(resultStSp_2)),
                  ylow = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,2], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,2],Y_test[,2]) )) + 
  ggtitle("AR1 model - 3rd order differences data Missouri, mixed ages, post high schooler")


Y_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],colMeans(resultStSp_1))
Y_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_2nddiff_pred_1)
Y_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_1stdiff_pred_1)

Y_low_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_low_2nddiff_pred_1)
Y_low_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_low_1stdiff_pred_1)

Y_up_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_up_2nddiff_pred_1)
Y_up_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_up_1stdiff_pred_1)

ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,1],Y_pred_1), 
                  ylow = c(Y[26,1],Y_low_pred_1), 
                  yup = c(Y[26,1],Y_up_pred_1)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Yh[1:26], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,1] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, high schooler")



Y_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],colMeans(resultStSp_2))
Y_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_2nddiff_pred_2)
Y_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_1stdiff_pred_2)

Y_low_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_low_2nddiff_pred_2)
Y_low_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_low_1stdiff_pred_2)

Y_up_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_up_2nddiff_pred_2)
Y_up_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_up_1stdiff_pred_2)

ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,2],Y_pred_2), 
                  ylow = c(Y[26,2],Y_low_pred_2), 
                  yup = c(Y[26,2],Y_up_pred_2)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y[1:26,2], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,2] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, post high schooler")

# Goodness of fit
llik   <- data.frame(rstan::extract(StSp_model2, "log_lik")[[1]])

p_WAIC_1 <- sum(apply(llik[,1:23], 2, var))
lppd_1   <- sum(apply(llik[,1:23], 2, function(x) log(mean(exp(x)))))
WAIC_1   <- - 2 * lppd_1 + 2 * p_WAIC_1 # 91.20442

p_WAIC_2 <- sum(apply(llik[,24:46], 2, var))
lppd_2   <- sum(apply(llik[,24:46], 2, function(x) log(mean(exp(x)))))
WAIC_2   <- - 2 * lppd_2 + 2 * p_WAIC_2 # 135.9375

# MSE
MSE_1 <- MSE(Y_pred_1, Y[27:29,1]) # 2.775141
MSE_2 <- MSE(Y_pred_2, Y[27:29,2]) # 2.121286

tab <- rbind(tab,
             c(WAIC_1, MSE_1, WAIC_2, MSE_2))
#-------------------------------------- END ----------------------------------#
###################################################################################


########################### THIRD MODEL: ARIMA(1,3,0) ############################
#------------------------------------------
# SHARED DATA  alpha_j|a0, s0  ~ N(a0, s0)   j = 1,2
#                phi_j|p0, g0  ~ N(p0, g0)   j = 1,2
# z_tj = alpha_j + phi_j * z_t-1,j + eps_j
data_StSp3 <-list(N = N_train,
                  A = A,
                  Y = Y_train,
                  alpha_m0 = 0,
                  alpha_phi = 0,
                  alpha_m0_var = 10,
                  alpha_phi_var = 10,
                  a_sigma2 = 2,
                  b_sigma2 = 10)

inits <- function() 
{
  list(m0 = c(0,0),
       phi = c(0,0), 
       sigma2 = c(var(Yh),var(Yp)),
       mu_m0 = 0,
       mu_phi = 0,
       sigma2m0 = 10,
       sigma2phi = 10)
}

# StSp_model3 <- rstan::stan(file = "StSp_model_2ages_third_version.stan",
#                    data = data_StSp3,
#                    chains = 2,
#                    iter = 30000,
#                    warmup = 10000,
#                    thin= 20,
#                    seed = 42,
#                    init = inits,
#                    algorithm = 'NUTS')
# 
# save(StSp_model3, file="AR1_mixed_ages_3rd_diff_3rd_version.dat")
load("AR1_mixed_ages_3rd_diff_3rd_version.dat")

rstan::traceplot(StSp_model3, pars = c("m0[1]", "m0[2]", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"), inc_warmup = TRUE)

# Diagnostic

# coda
coda_chain <- rstan::As.mcmc.list(StSp_model3, pars = c("m0[1]", "m0[2]", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"))
summary(coda_chain)
#              Mean     SD Naive SE Time-series SE
# m0[1]      0.0542 0.3940 0.008810       0.008651       good
# m0[2]      0.4846 0.8800 0.019678       0.019145       good
# phi[1]    -0.8409 0.1383 0.003094       0.003094       good
# phi[2]    -0.5752 0.1890 0.004227       0.004227       good
# sigma2[1]  3.4437 1.0561 0.023614       0.022808       good
# sigma2[2] 17.1680 5.5467 0.124028       0.114457       good

#              2.5%      25%      50%     75%   97.5%
# m0[1]     -0.7303 -0.21145  0.05681  0.3057  0.8146
# m0[2]     -1.2247 -0.09782  0.47559  1.0456  2.2617
# phi[1]    -1.1257 -0.93271 -0.83736 -0.7480 -0.5833    NOT true that |phi| < 1 in a 95% CI
# phi[2]    -0.9610 -0.69797 -0.57312 -0.4544 -0.1937    true that |phi| < 1 in a 95% CI
# sigma2[1]  1.9440  2.69624  3.27059  3.9724  6.0468
# sigma2[2]  9.5735 13.29514 16.05867 19.9049 30.6383

# Gelman and rubin's convergence diagnostic
coda::gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
#ok, parameters ~ Normal have 1 as an estimate

# Geweke's convergence diagnostic
coda::geweke.diag(coda_chain, frac1=0.1, frac2=0.5)

# Autocorrelation and plot
rstan::stan_ac(StSp_model3, pars = c("m0[1]", "m0[2]", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"))
# ok, indepentent samples

# Posterior distributions
plot_post <- StSp_model3 %>% 
  rstan::extract(c("m0[1]", "m0[2]", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")


# Predictions
params <- data.frame(rstan::extract(StSp_model3, c("m0[1]", "m0[2]", "phi[1]", "phi[2]", "sigma2[1]", "sigma2[2]"), perm = T))
resultStSp_1 <- matrix(NA, nrow = nrow(params), ncol = N_test)
resultStSp_2 <- matrix(NA, nrow = nrow(params), ncol = N_test)
for (i in 1:N_test) {
  if(i == 1){
    resultStSp_1[,i] <- apply(params, 1, function(x) rnorm(1, x[1] + x[3] * Y_train[N_train,1], sqrt(x[5])))
    resultStSp_2[,i] <- apply(params, 1, function(x) rnorm(1, x[2] + x[4] * Y_train[N_train,2], sqrt(x[6])))
  } else {
    resultStSp_1[,i] <- apply(cbind(params, resultStSp_1[,i-1]), 1, function(x) rnorm(1, x[1] + x[3] * x[7], sqrt(x[5])))
    resultStSp_2[,i] <- apply(cbind(params, resultStSp_2[,i-1]), 1, function(x) rnorm(1, x[2] + x[4] * x[7], sqrt(x[6])))
  }
}

ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,1],colMeans(resultStSp_1)),
                  ylow = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,1], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,1],Y_test[,1]) )) + 
  ggtitle("AR1 model - 3rd order differences data Missouri, mixed ages, high schooler")

ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,2],colMeans(resultStSp_2)),
                  ylow = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,2], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,2],Y_test[,2]) )) + 
  ggtitle("AR1 model - 3rd order differences data Missouri, mixed ages, post high schooler")


Y_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],colMeans(resultStSp_1))
Y_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_2nddiff_pred_1)
Y_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_1stdiff_pred_1)

Y_low_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_low_2nddiff_pred_1)
Y_low_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_low_1stdiff_pred_1)

Y_up_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_up_2nddiff_pred_1)
Y_up_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_up_1stdiff_pred_1)

ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,1],Y_pred_1), 
                  ylow = c(Y[26,1],Y_low_pred_1), 
                  yup = c(Y[26,1],Y_up_pred_1)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Yh[1:26], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,1] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, high schooler")



Y_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],colMeans(resultStSp_2))
Y_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_2nddiff_pred_2)
Y_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_1stdiff_pred_2)

Y_low_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_low_2nddiff_pred_2)
Y_low_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_low_1stdiff_pred_2)

Y_up_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_up_2nddiff_pred_2)
Y_up_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_up_1stdiff_pred_2)

ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,2],Y_pred_2), 
                  ylow = c(Y[26,2],Y_low_pred_2), 
                  yup = c(Y[26,2],Y_up_pred_2)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y[1:26,2], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,2] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, post high schooler")

# Goodness of fit
llik   <- data.frame(rstan::extract(StSp_model3, "log_lik")[[1]])

p_WAIC_1 <- sum(apply(llik[,1:23], 2, var))
lppd_1   <- sum(apply(llik[,1:23], 2, function(x) log(mean(exp(x)))))
WAIC_1   <- - 2 * lppd_1 + 2 * p_WAIC_1 #91.59127

p_WAIC_2 <- sum(apply(llik[,24:46], 2, var))
lppd_2   <- sum(apply(llik[,24:46], 2, function(x) log(mean(exp(x)))))
WAIC_2   <- - 2 * lppd_2 + 2 * p_WAIC_2 # 136.1916


# MSE for evaluate forcasting performance
MSE_1 <- MSE(Y_pred_1, Y[27:29,1]) # 1.990304
MSE_2 <- MSE(Y_pred_2, Y[27:29,2]) # 5.683792
tab <- rbind(tab,
             c(WAIC_1, MSE_1, WAIC_2, MSE_2))
#--------------------------------------- END ----------------------------------#
##################################################################################


########################### FOURTH MODEL: ARIMA(2,3,3) ############################
#------------------------------------------
# d_tj = alpha + phi1_j * d_t-1,j + phi2_j * d_t-2,j 
#              + beta1_j * eps_t-1,j +beta2_j * eps_t-2,j+ beta3_j * eps_t-3,j
#              + eps_t


data_StSp <-list(N = N_train,
                 A = A,
                 Y = Y_train,
                 a_sigma2 = 2,
                 b_sigma2 = 10,
                 sigma2m0 = 10)

inits <- function() 
{
  list(phi1 = c(0,0),
       phi2 = c(0,0),
       beta1 = c(0,0),
       beta2 = c(0,0),
       beta3 = c(0,0),
       m0 = 0, 
       sigma2 = c(var(Yh),var(Yp)))
}

# StSp_model <- rstan::stan(file = "ARIMA_2_3_3.stan",
#                           data = data_StSp,
#                           chains = 2,
#                           iter = 30000,
#                           warmup = 10000,
#                           thin= 20,
#                           seed = 42,
#                           init = inits,
#                           algorithm = 'NUTS')
# 
# save(StSp_model, file="ARIMA_2_3_3_mixed_ages_3rd_diff.dat")
load("ARIMA_2_3_3_mixed_ages_3rd_diff.dat")

rstan::traceplot(StSp_model, pars = c("m0", "phi1[1]", "phi1[2]", "phi2[1]", "phi2[2]", "beta1[1]", "beta1[2]", "beta2[1]", "beta2[2]", "beta3[1]", "beta3[2]", "sigma2[1]", "sigma2[2]"), inc_warmup = TRUE)

# Diagnostic

# coda
coda_chain <- rstan::As.mcmc.list(StSp_model, pars = c("m0", "phi1[1]", "phi1[2]", "phi2[1]", "phi2[2]", "beta1[1]", "beta1[2]", "beta2[1]", "beta2[2]", "beta3[1]", "beta3[2]", "sigma2[1]", "sigma2[2]"))
summary(coda_chain) # quite ugly

# Gelman and rubin's convergence diagnostic
coda::gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
coda::geweke.diag(coda_chain, frac1=0.1, frac2=0.5)

# Autocorrelation and plot
rstan::stan_ac(StSp_model, pars = c("m0", "phi1[1]", "phi1[2]", "phi2[1]", "phi2[2]", "beta1[1]", "beta1[2]", "beta2[1]", "beta2[2]", "beta3[1]", "beta3[2]", "sigma2[1]", "sigma2[2]"))
# too much correlation even if we have thinned the chain!

# Posterior distributions
plot_post <- StSp_model %>% 
  rstan::extract(c("m0", "phi1[1]", "phi1[2]", "phi2[1]", "phi2[2]", "beta1[1]", "beta1[2]", "beta2[1]", "beta2[2]", "beta3[1]", "beta3[2]", "sigma2[1]", "sigma2[2]")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# Predictions
# m0 + phi1[a] * Y[t-1,a] + phi2[a] * Y[t-2,a]  + beta1[a] * epsilon[t-1,a] + beta2[a] * epsilon[t-2,a]+beta3[a] * epsilon[t-3,a]
epsilon   <- data.frame(rstan::extract(StSp_model, "epsilon"))

params <- data.frame(rstan::extract(StSp_model, c("m0", "phi1[1]", "phi1[2]", "phi2[1]", "phi2[2]", "beta1[1]", "beta1[2]", "beta2[1]", "beta2[2]", "beta3[1]", "beta3[2]", "sigma2[1]", "sigma2[2]"), perm = T))

resultStSp_1 <- matrix(NA, nrow = nrow(params), ncol = N_test)
resultStSp_2 <- matrix(NA, nrow = nrow(params), ncol = N_test)


resultStSp_1[,1] <- apply(cbind(params,epsilon), 1, function(x) rnorm(1, x[1] + x["phi1.1."] * Y_train[N_train,1]+x["phi2.1."] * Y_train[N_train-1,1] + x["beta1.1."] * x["epsilon.23.1"] + x["beta2.1."] * x["epsilon.22.1"] + x["beta3.1."] * x["epsilon.21.1"] , sqrt(x["sigma2.1."])))
resultStSp_2[,1] <- apply(cbind(params,epsilon), 1, function(x) rnorm(1, x[1] + x["phi1.2."] * Y_train[N_train,2]+x["phi2.2."] * Y_train[N_train-1,2] + x["beta1.2."] * x["epsilon.23.2"] + x["beta2.2."] * x["epsilon.22.2"] + x["beta3.2."] * x["epsilon.21.2"] , sqrt(x["sigma2.2."])))

# for 2nd year prediction, we don't have true value for 1st prediction.so we cannot compute the eplison_t-1 and I just set it as 0;
resultStSp_1[,2] <- apply(cbind(params,epsilon, resultStSp_1[,1]), 1, function(x) rnorm(1, x[1] + x["phi1.1."] * x[60]+x["phi2.1."] * Y_train[N_train,1] + x["beta1.1."] * 0  + x["beta2.1."] * x["epsilon.23.1"] + x["beta3.1."] * x["epsilon.22.1"] , sqrt(x["sigma2.1."])))
resultStSp_2[,2] <- apply(cbind(params,epsilon, resultStSp_2[,1]), 1, function(x) rnorm(1, x[1] + x["phi1.2."] * x[60]+x["phi2.2."] * Y_train[N_train,2] + x["beta1.2."] * 0  + x["beta2.2."] * x["epsilon.23.2"] + x["beta3.2."] * x["epsilon.22.2"] , sqrt(x["sigma2.2."])))

resultStSp_1[,3] <- apply(cbind(params,epsilon, resultStSp_1[,1], resultStSp_1[,2]), 1, function(x) rnorm(1, x[1] + x["phi1.1."] * x[61]+x["phi2.1."] * x[60] + x["beta1.1."] * 0  + x["beta2.1."] * 0 + x["beta3.1."] * x["epsilon.23.1"] , sqrt(x["sigma2.1."])))
resultStSp_2[,3] <- apply(cbind(params,epsilon, resultStSp_2[,1], resultStSp_1[,2]), 1, function(x) rnorm(1, x[1] + x["phi1.2."] * x[61]+x["phi2.2."] * x[60] + x["beta1.2."] * 0  + x["beta2.2."] * 0 + x["beta3.2."] * x["epsilon.23.2"] , sqrt(x["sigma2.2."])))


ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,1],colMeans(resultStSp_1)),
                  ylow = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,1],apply(resultStSp_1, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,1], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,1],Y_test[,1]) )) + 
  ggtitle("AR1 model - Third order differences data Missouri, mixed ages, high schooler")



ggplot(data.frame(x = (N_train) :(N_train + N_test), 
                  y = c(Y_train[N_train,2],colMeans(resultStSp_2)),
                  ylow = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.05)),
                  yup = c(Y_train[N_train,2],apply(resultStSp_2, 2, quantile, p = 0.95)))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y_train[,2], x = 1:N_train), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = N_train), col = 2, lty = 2) +
  geom_point(aes(x=(N_train) :(N_train + N_test), y=c(Y_train[N_train,2],Y_test[,2]) )) + 
  ggtitle("AR1 model - Third order differences data Missouri, mixed ages, post high schooler")


Y_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],colMeans(resultStSp_1))
Y_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_2nddiff_pred_1)
Y_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_1stdiff_pred_1)

Y_low_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_low_2nddiff_pred_1)
Y_low_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_low_1stdiff_pred_1)

Y_up_2nddiff_pred_1<-diff_to_timeseries(Y_2nddiff[N_train+1,1],apply(resultStSp_1, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_1<-diff_to_timeseries(Y_1stdiff[N_train+2,1],Y_up_2nddiff_pred_1)
Y_up_pred_1        <-diff_to_timeseries(Y[N_train+3,1],Y_up_1stdiff_pred_1)

ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,1],Y_pred_1), 
                  ylow = c(Y[26,1],Y_low_pred_1), 
                  yup = c(Y[26,1],Y_up_pred_1)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Yh[1:26], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,1] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, high schooler")



Y_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],colMeans(resultStSp_2))
Y_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_2nddiff_pred_2)
Y_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_1stdiff_pred_2)

Y_low_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.05))
Y_low_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_low_2nddiff_pred_2)
Y_low_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_low_1stdiff_pred_2)

Y_up_2nddiff_pred_2<-diff_to_timeseries(Y_2nddiff[N_train+1,2],apply(resultStSp_2, 2, quantile, p = 0.95))
Y_up_1stdiff_pred_2<-diff_to_timeseries(Y_1stdiff[N_train+2,2],Y_up_2nddiff_pred_2)
Y_up_pred_2        <-diff_to_timeseries(Y[N_train+3,2],Y_up_1stdiff_pred_2)

ggplot(data.frame(x = 26:29, 
                  y = c(Y[26,2],Y_pred_2), 
                  ylow = c(Y[26,2],Y_low_pred_2), 
                  yup = c(Y[26,2],Y_up_pred_2)))+
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y[1:26,2], x = 1:26), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = 26), col = 2, lty = 2) +
  geom_point(aes(x=c(26,27,28,29), y=Y[26:29,2] )) + 
  ggtitle("AR1 model - BIRTH rate data Missouri, mixed ages, post high schooler")


# Goodness of fit
llik   <- data.frame(rstan::extract(StSp_model, "log_lik")[[1]])

p_WAIC_1 <- sum(apply(llik[,4:23], 2, var))
lppd_1   <- sum(apply(llik[,4:23], 2, function(x) log(mean(exp(x)))))
WAIC_1   <- - 2 * lppd_1 + 2 * p_WAIC_1   ##73.4827

p_WAIC_2 <- sum(apply(llik[,27:46], 2, var))
lppd_2   <- sum(apply(llik[,27:46], 2, function(x) log(mean(exp(x)))))
WAIC_2   <- - 2 * lppd_2 + 2 * p_WAIC_2   ##104.6212


# MSE for evaluate forcasting performance
MSE_1 <- MSE(Y_pred_1, Y[27:29,1])  #0.1279609
MSE_2 <- MSE(Y_pred_2, Y[27:29,2])  #2.797991
tab <- rbind(tab,
             c(WAIC_1, MSE_1, 0.03884597, WAIC_2, MSE_2, 0.03202568))

#----------------------------------------- END ----------------------------------#
###################################################################################


TAB <- data.frame(tab)
colnames(TAB) <- c("WAIC High", "MSE High", "WAIC PostHigh", "MSE PostHigh")
row.names(TAB) <- c("ARIMA(1,3,0) 1st version", "ARIMA(1,3,0) 2nd version", "ARIMA(1,3,0) 3rd version", "ARIMA(2,3,3)")
TAB


#---------------------------- TEST ----------------------------
#              Mean     SD Naive SE Time-series SE
# m0         0.1228 0.3502 0.007830       0.007827      good
# phi[1]    -0.8385 0.1377 0.003080       0.003122      good
# phi[2]    -0.5716 0.1828 0.004087       0.004058      good
# sigma2[1]  3.3809 1.0522 0.023527       0.025329      ok
# sigma2[2] 16.6686 5.1460 0.115067       0.109909      ok
set.seed(1441)
state <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", 
           "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia", 
           "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky",
           "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
           "Mississippi", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey",
           "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma",
           "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
           "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia",
           "Wisconsin", "Wyoming")
N_try = 25
p <- c(0.1228, -0.8385, -0.5716, 3.3809, 16.6686)

### run directly load(file_names.dat) for the result
# mse_tot <- matrix(NA, nrow = length(state), ncol = 2)
# ts_high_all <- matrix(NA, nrow = length(state), ncol = 29)
# ts_post_all <- matrix(NA, nrow = length(state), ncol = 29)
# ts_high_diff <- matrix(NA, nrow = length(state), ncol = 26)
# ts_post_diff <- matrix(NA, nrow = length(state), ncol = 26)
# ts_high_pred <- matrix(NA, nrow = length(state), ncol = 26)
# ts_post_pred <- matrix(NA, nrow = length(state), ncol = 26)
 

load("mse_usa.dat")
load("ts_high_all.dat")
load("ts_post_all.dat")
load("ts_high_diff.dat")
load("ts_post_diff.dat")
load("ts_high_pred.dat")
load("ts_post_diff.dat")

# for (s in 1:length(state)) {
#   ts_high_all[s,] <- Yh_s <- BIRTH_data$State.Rate[BIRTH_data$State == state[s] & BIRTH_data$Age.Group..Years. == "15-17 years"]
#   ts_post_all[s,] <- Yp_s <- BIRTH_data$State.Rate[BIRTH_data$State == state[s] & BIRTH_data$Age.Group..Years. == "18-19 years"]
#   Y_real <- cbind(Yh_s, Yp_s)
#   
#   Y_3rddiff <- diff(Y_real,diff=3)
#   ts_high_diff[s,] <- Y_3rddiff[,1]
#   ts_post_diff[s,] <- Y_3rddiff[,2]
#   Y_out_s <- Y_3rddiff[1,]
#   res_h <- rep(NA, N_try)
#   res_p <- rep(NA, N_try)
#   
#   for (i in 1:N_try) {
#     if(i == 1){
#       res_h[i] <- rnorm(1, p[1] + p[2] * Y_out_s[1], sqrt(p[4]))
#       res_p[i] <- rnorm(1, p[1] + p[3] * Y_out_s[2], sqrt(p[5]))
#     } else {
#       res_h[i] <- rnorm(1, p[1] + p[2] * res_h[i-1], sqrt(p[4]))
#       res_p[i] <- rnorm(1, p[1] + p[3] * res_p[i-1], sqrt(p[5]))
#     }
#   }
#   
#   Y_out_s <- rbind(Y_out_s, cbind(res_h,res_p))
#   ts_high_pred[s,] <- Y_out_s[,1]
#   ts_post_pred[s,] <- Y_out_s[,2]
#   mse_h <- MSE(Y_out_s[,1], Y_3rddiff[,1])
#   mse_p <- MSE(Y_out_s[,2], Y_3rddiff[,2])
#   mse_tot[s,] <- c(mse_h, mse_p)
#   
# }
# 
# row.names(mse_tot) <- state
# colnames(mse_tot) <- c("mse high", "mse post-high")
# save(mse_tot, file="mse_usa.dat")
# save(ts_high_all, file="ts_high_all.dat")
# save(ts_post_all, file="ts_post_all.dat")
# save(ts_high_diff, file="ts_high_diff.dat")
# save(ts_post_diff, file="ts_post_diff.dat")
# save(ts_high_pred, file="ts_high_pred.dat")
# save(ts_post_pred, file="ts_post_pred.dat")

mse_tot

m_high <- min(mse_tot[,1]) #Virginia
m_post <- min(mse_tot[,2]) #Ohio 


# Minimum mse:
ggplot(data.frame(x = 1:26), y = ts_high_pred[mse_tot[,1] == m_high,]) +
  geom_line(data = data.frame(Y = ts_high_pred[mse_tot[,1] == m_high,], x = 1:26), 
            aes(x = x, y = Y), lty = 2, col = "red") +
  geom_line(data = data.frame(Y = ts_high_diff[mse_tot[,1] == m_high,], x = 1:26), 
            aes(x = x, y = Y), col = "blue") + theme_bw() +
  ggtitle("Forecasting on State Virginia for Highschool")


ggplot(data.frame(x = 1:26), y = ts_post_diff[mse_tot[,2] == m_post,]) +
  geom_line(data = data.frame(Y = ts_post_diff[mse_tot[,2] == m_post,], x = 1:26), 
            aes(x = x, y = Y), lty = 2, col = "red") +
  geom_line(data = data.frame(Y = ts_post_pred[mse_tot[,2] == m_post,], x = 1:26), 
            aes(x = x, y = Y), col = "blue") + theme_bw() +
  ggtitle("Forecasting on State Ohio for Post-Highschool")


# BAD mse but state is close to missouri
ggplot(data.frame(x = 1:26), y = ts_post_diff[state == "Nebraska",]) +
  geom_line(data = data.frame(Y = ts_post_diff[state == "Nebraska",], x = 1:26), 
            aes(x = x, y = Y), lty = 2, col = "red") +
  geom_line(data = data.frame(Y = ts_post_pred[state == "Nebraska",], x = 1:26), 
            aes(x = x, y = Y), col = "blue") + theme_bw() +
  ggtitle("Forecasting on State Nebraska for Post-Highschool")


# GOOD mse but state is far away from missouri
ggplot(data.frame(x = 1:26), y = ts_high_pred[state == "Oregon",]) +
  geom_line(data = data.frame(Y = ts_high_pred[state == "Oregon",], x = 1:26), 
            aes(x = x, y = Y), lty = 2, col = "red") +
  geom_line(data = data.frame(Y = ts_high_diff[state == "Oregon",], x = 1:26), 
            aes(x = x, y = Y), col = "blue") + theme_bw() +
  ggtitle("Forecasting on State Oregon for Highschool")

# the most weird mse
ggplot(data.frame(x = 1:26), y = ts_high_pred[state == "District of Columbia",]) +
  geom_line(data = data.frame(Y = ts_high_pred[state == "District of Columbia",], x = 1:26), 
            aes(x = x, y = Y), lty = 2, col = "red") +
  geom_line(data = data.frame(Y = ts_high_diff[state == "District of Columbia",], x = 1:26), 
            aes(x = x, y = Y), col = "blue") + theme_bw() +
  ggtitle("Forecasting on State District of Columbia for Highschool")

graphics.off()
rm(list=ls())

