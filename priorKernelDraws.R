library(rstan)
library(cmdstanr)
library(bayesplot)
library(plyr)
library(posterior)
library(fields)
library(reticulate)
library(MASS)
library(tidyverse)
library(patchwork)
library(TruncatedNormal)
#-----------------------------------------
#Setting params
#-----------------------------------------

sigma_est = 1

length_scale1_est = 4

length_scale2_est = 1

length_scale3_est = 0.5

length_scale4_est = 1.7

period_est = 0.75

input1 = seq(0,10,length.out = 1000)

#-----------------------------------------
#Quasi-periodic kernel draws
#-----------------------------------------

Cov_quasi_period <- function(x1, x1p, sigma = sigma_est, length_scale1 = length_scale1_est,
                       length_scale2 = length_scale2_est, period = period_est){
  res = sigma^2 * 
    (exp(-((outer(x1, x1p, "-"))^2)/(2*(length_scale1)^2)) *
      exp(-2/(length_scale2)^2 * (sin(pi * abs(outer(x1, x1p, "-"))/period))^2))
  return(res)
}

prior_cov_matrix_QP = Cov_quasi_period(input1,input1)

prior_sample_QP = mvrnorm(n = 5, mu = rep(0,length(input1)), Sigma = prior_cov_matrix_QP)

#-----------------------------------------
#Squared exponential kernel draws
#-----------------------------------------

Cov_SE <- function(x1, x1p, sigma = sigma_est, length_scale3 = length_scale3_est){
  res = sigma^2 * 
    exp(-((outer(x1, x1p, "-"))^2)/(2*(length_scale3)^2))
  return(res)
}

prior_cov_matrix_SE = Cov_SE(input1,input1)

prior_sample_SE = mvrnorm(n = 5, mu = rep(0,length(input1)), Sigma = prior_cov_matrix_SE)

#-----------------------------------------
#Periodic kernel draws
#-----------------------------------------

Cov_periodic <- function(x1, x1p, sigma = sigma_est,length_scale4 = length_scale4_est,
                             period = period_est){
  res = sigma^2 * 
       (exp(-2/(length_scale4)^2 * (sin(pi * abs(outer(x1, x1p, "-"))/period))^2))
  return(res)
}

prior_cov_matrix_periodic = Cov_periodic(input1,input1)

prior_sample_periodic = mvrnorm(n = 5, mu = rep(0,length(input1)), Sigma = prior_cov_matrix_periodic)

#-----------------------------------------
#Plotting
#-----------------------------------------

p1 = ggplot() + 
  geom_line(aes(x = input1, y = prior_sample_periodic[1,]), col = "#022020", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_periodic[2,] + 3), col = "#003c3c", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_periodic[3,] + 6), col = "#007d69", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_periodic[4,] + 9), col = "#00a87e", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_periodic[5,] + 12), col = "#00c896", lwd = 0.8)+
  geom_vline(xintercept = seq(0,10,.75), alpha = .2) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title = "Periodic kernel",
       y = "",
       x = "")

p2 = ggplot() + 
  geom_line(aes(x = input1, y = prior_sample_SE[1,]), col = "red", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_SE[2,]), col = "black", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_SE[3,]), col = "blue", lwd = 0.8)+
  labs(title = "Squared exponential kernel draws",
       y = "Output",
       x = "Input")

p2

p3 = ggplot() + 
  geom_line(aes(x = input1, y = prior_sample_QP[1,]), col = "#022020", lwd = 0.8) +
  geom_line(aes(x = input1, y = prior_sample_QP[2,] + 3), col = "#003c3c", lwd = 0.8) +
  geom_line(aes(x = input1, y = prior_sample_QP[3,] + 6), col = "#007d69", lwd = 0.8) +
  geom_line(aes(x = input1, y = prior_sample_QP[4,] + 9), col = "#00a87e", lwd = 0.8) +
  geom_line(aes(x = input1, y = prior_sample_QP[5,] + 12), col = "#00c896", lwd = 0.8) +
  geom_vline(xintercept = seq(0,10,.75), alpha = .2) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title = "Quasi-periodic kernel (squared exponential x periodic)",
       y = "",
       x = "Input")

p1/p3

p2

observed_x = c(1,7.75)
observed_y = c(3,-1.5)

pred_mean <- function(x1_pred, x1_obs, y, nugget = 0){
  C_xX <- Cov_SE(x1_pred, x1_obs)
  C_XX <- Cov_SE(x1_obs, x1_obs)
  m <- C_xX %*% solve(C_XX + (nugget*diag(ncol(C_XX)))) %*% y
  return(m)
}

pred_cov <- function(x1_pred, x1_obs, nugget = 0){
  C_xX <- Cov_SE(x1_pred, x1_obs)
  C_XX <- Cov_SE(x1_obs, x1_obs)
  s <- Cov_SE(x1_pred, x1_pred) - C_xX %*% solve(C_XX + (nugget*diag(ncol(C_XX)))) %*% t(C_xX)
  return(s)
}

mx = pred_mean(input1, observed_x, observed_y)
s2x = pred_cov(input1, observed_x) 

conditional_draws = mvrnorm(n = 3, mu = mx, Sigma = s2x)

ggplot() + 
  geom_point(aes(x = observed_x, y = observed_y), col = "#007d69", cex = 3) + 
  xlim(0,10) + 
  ylim(-4,4)

ggplot() + 
  geom_line(aes(x = input1, conditional_draws[1,]), col = "red", lwd = 0.8) + 
  geom_line(aes(x = input1, conditional_draws[2,]), col = "black", lwd = 0.8) + 
  geom_line(aes(x = input1, conditional_draws[3,]), col = "blue", lwd = 0.8) + 
  geom_point(aes(x = observed_x, y = observed_y), col = "#007d69", cex = 3) + 
  xlim(0,10) + 
  ylim(-4,4)
