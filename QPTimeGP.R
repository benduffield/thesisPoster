library(rstan)
library(cmdstanr)
library(bayesplot)
library(plyr)
library(posterior)
library(fields)
library(reticulate)
library(MASS)

#-----------------------------------------
#Creating model
#-----------------------------------------

#import .stan file
playground_model = cmdstan_model(stan_file = "GP_periodic_time.stan")

#run CVS model for 0.75s with 30 time steps
Vspt_obs = Vspt
Vrv_obs = Vrv
Vlv_obs = Vlv
tsteps_obs = tsteps

mydata = list(
  N = length(Vrv_obs),
  x1 = Vrv_obs,
  x2 = Vlv_obs,
  x3 = tsteps_obs,
  y = Vspt_obs
)

#-----------------------------------------
#MAP estimates for parameters (using 30 points generated by NR in 0-0.75 second range)
#-----------------------------------------

MAP_params = playground_model$optimize(data = mydata, jacobian = TRUE)

MAP_params$summary()

sigma_MAP = as.double(MAP_params$summary("sigma")[2])

length_scale1_MAP = as.double(MAP_params$summary("length_scale1")[2])

length_scale2_MAP = as.double(MAP_params$summary("length_scale2")[2])

length_scale3_MAP = as.double(MAP_params$summary("length_scale3")[2])

period_MAP = 0.75

tau_sq_MAP = 1e-7

#Posterior sampling and plots
param_samples = playground_model$sample(data = mydata,
                                        chains = 3,
                                        parallel_chains = 3,
                                        iter_warmup = 1000,
                                        iter_sampling = 1000)

mcmc_trace(param_samples$draws(),pars = c("sigma", "length_scale1", "length_scale2","length_scale3"))
param_samples
param_samples$draws()
mcmc_hist(param_samples$draws())

#-----------------------------------------
#Setting up mean and covariance for prediction
#-----------------------------------------

#QP kernel function
Cov_QP <- function(x1, x1p, x2, x2p, x3, x3p, sigma = sigma_MAP, length_scale1 = length_scale1_MAP,
                   length_scale2 = length_scale2_MAP, length_scale3 = length_scale3_MAP,
                   period = period_MAP){
  
  N1 = length(x1)
  N2 = length(x1p)
  
  K = matrix(, nrow = N1, ncol = N2)
  
  for (i in 1:N1) {
    for (j in 1:N2) {
      K[i, j] = ((exp(-(x1[i] - x1p[j])^2 / (2 * length_scale1^2))) *
                   (exp(-(x2[i] - x2p[j])^2 / (2 * length_scale2^2))) * 
                   (exp((-2/(length_scale3)^2) * (sin(pi * abs(x3[i]-x3p[j])/period))^2))) *
        sigma^2
    }
  }
  
  return(K)
}


#Posterior mean function
pred_mean <- function(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs, y, nugget = tau_sq_MAP){
  C_pred_obs <- Cov_QP(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs)
  C_obs_obs <- Cov_QP(x1_obs, x1_obs, x2_obs, x2_obs, x3_obs, x3_obs)
  m <-C_pred_obs %*% solve(C_obs_obs + (nugget*diag(ncol(C_obs_obs)))) %*% y
  return(m)
}

#Posterior covariance matrix
pred_cov_QP <- function(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs, nugget = tau_sq_MAP){
  C_pred_obs <- Cov_QP(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs)
  C_obs_obs <- Cov_QP(x1_obs, x1_obs, x2_obs, x2_obs, x3_obs, x3_obs)
  C_pred_pred = Cov_QP(x1_pred, x1_pred, x2_pred, x2_pred, x3_pred, x3_pred)
  s <- C_pred_pred - (C_pred_obs %*% solve(C_obs_obs + (nugget*diag(ncol(C_obs_obs)))) %*% t(C_pred_obs))
  print(isSymmetric(C_pred_pred))
  print(isSymmetric((C_pred_obs %*% solve(C_obs_obs + (nugget*diag(ncol(C_obs_obs)))) %*% t(C_pred_obs))))
  print((C_pred_obs %*% solve(C_obs_obs + (nugget*diag(ncol(C_obs_obs)))) %*% t(C_pred_obs))[1,5])
  print((C_pred_obs %*% solve(C_obs_obs + (nugget*diag(ncol(C_obs_obs)))) %*% t(C_pred_obs))[5,1])
  return(s)
}

#-----------------------------------------
#Predicted values (evaluate NR at 1000 points in 0-10 second range)
#-----------------------------------------

#Run CVS model for longer with more timesteps
Vrv_pred = Vrv
Vlv_pred = Vlv
time_pred = tsteps

#DF where column 1 = mean, 2 = var, 3 = +3 s.d., 4 = -3 s.d.
Pred_df = matrix(nrow = length(Vrv_pred), ncol = 4)

Pred_df[,1] = pred_mean(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs, Vspt_obs)
Pred_df[,2] = diag(pred_cov_QP(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs))
Pred_df[,3] = Pred_df[,1] + (3 * sqrt(Pred_df[,2]))
Pred_df[,4] = Pred_df[,1] - (3 * sqrt(Pred_df[,2]))

#Plotting
QP_plot = ggplot() + 
  geom_line(aes(x = time_pred, y = Pred_df[,1]), col = "#007d69", lwd = 0.8) + 
  geom_point(aes(x = tsteps_obs, y = Vspt_obs), col = "#003c3c") +
  geom_line(aes(x = time_pred, y = Pred_df[,3]), lty = 2, col = "#813534") + 
  geom_line(aes(x = time_pred, y = Pred_df[,4]), lty = 2, col = "#813534") +
  labs(title = "QP kernel interpolation",
       x = "Time",
       y = "Vspt")

QP_plot

#Compute mean vector and cov matrix for MVN
Cov_matrix = pred_cov_QP(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs)
Mean_vector = pred_mean(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs, Vspt_obs)

View(Cov_matrix) #NOT SYMMETRIC

#sample
posterior_samples = TruncatedNormal::rtmvnorm(n = 1, mu = Mean_vector, sigma = Cov_matrix,
                                              lb = rep(0, length(Mean_vector)))


#Plot samples
p1 = ggplot() + 
  geom_line(aes(x = time_pred, y = posterior_samples[1,]), col = "#022020", lwd = 0.8)+
  geom_line(aes(x = time_pred, y = posterior_samples[2,] + 3), col = "#003c3c", lwd = 0.8)+
  geom_line(aes(x = time_pred, y = posterior_samples[3,] + 6), col = "#007d69", lwd = 0.8)+
  geom_line(aes(x = time_pred, y = posterior_samples[4,] + 9), col = "#00a87e", lwd = 0.8)+
  geom_line(aes(x = time_pred, y = posterior_samples[5,] + 12), col = "#00c896", lwd = 0.8)+
  geom_vline(xintercept = seq(0,10,.75), alpha = .2) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title = "Periodic kernel",
       y = "",
       x = "")
p1  

eigen(Cov_matrix, symmetric = TRUE)$values
