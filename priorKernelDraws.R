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
#-----------------------------------------
#Setting params
#-----------------------------------------

sigma_est = 1

length_scale1_est = 4

length_scale2_est = 1

length_scale3_est = 0.1

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
  geom_line(aes(x = input1, y = prior_sample_SE[1,]), col = "#022020", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_SE[2,] + 3), col = "#003c3c", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_SE[3,] + 6), col = "#007d69", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_SE[4,] + 9), col = "#00a87e", lwd = 0.8)+
  geom_line(aes(x = input1, y = prior_sample_SE[5,] + 12), col = "#00c896", lwd = 0.8)+
  geom_vline(xintercept = seq(0,10,.75), alpha = .2) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title = "Squared exponential kernel",
       y = "Output",
       x = "")

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

p1 / p2 / p3
  
