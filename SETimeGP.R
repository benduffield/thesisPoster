library(rstan)
library(cmdstanr)
library(bayesplot)
library(plyr)
library(posterior)
library(fields)
library(reticulate)
#-----------------------------------------
#Creating model
#-----------------------------------------

playground_model = cmdstan_model(stan_file = "GP_time.stan")

Vspt_obs = Vspt
Vrv_obs = Vrv
Vlv_obs = Vlv
tsteps_obs = tsteps

mydata = list(
  N = length(Vspt_obs),
  x1 = Vrv_obs,
  x2 = Vlv_obs,
  x3 = tsteps_obs,
  y = Vspt_obs
)

#-----------------------------------------
#MAP estimates for parameters
#-----------------------------------------

MAP_params = playground_model$optimize(data = mydata, jacobian = TRUE)

MAP_params$summary()

sigma_MAP_SE = as.double(MAP_params$summary("sigma")[2])

length_scale1_MAP_SE = as.double(MAP_params$summary("length_scale1")[2])

length_scale2_MAP_SE = as.double(MAP_params$summary("length_scale2")[2])

length_scale3_MAP_SE = as.double(MAP_params$summary("length_scale3")[2])

tau_sq_MAP_SE = 1e-7

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

Cov_Exp_sq <- function(x1, x1p, x2, x2p, x3, x3p, sigma = sigma_MAP_SE, length_scale1 = length_scale1_MAP_SE,
                       length_scale2 = length_scale2_MAP_SE, length_scale3 = length_scale3_MAP_SE){
  res = sigma^2 * 
    exp(-((outer(x1, x1p, "-"))^2)/(2*(length_scale1)^2)) * 
    exp(-((outer(x2, x2p, "-"))^2)/(2*(length_scale2)^2)) *
    exp(-((outer(x3, x3p, "-"))^2)/(2*(length_scale3)^2))
  return(res)
}


pred_mean_SE <- function(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs, y, nugget = tau_sq_MAP){
  C_xX <- Cov_Exp_sq(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs)
  C_XX <- Cov_Exp_sq(x1_obs, x1_obs, x2_obs, x2_obs, x3_obs, x3_obs)
  m <-C_xX %*% solve(C_XX + (nugget*diag(ncol(C_XX)))) %*% y
  return(m)
}

pred_cov_SE <- function(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs, nugget = tau_sq_MAP){
  C_xX <- Cov_Exp_sq(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs)
  C_XX <- Cov_Exp_sq(x1_obs, x1_obs, x2_obs, x2_obs, x3_obs, x3_obs)
  s <- Cov_Exp_sq(x1_pred, x1_pred, x2_pred, x2_pred, x3_pred, x3_pred) - C_xX %*% solve(C_XX + (nugget*diag(ncol(C_XX)))) %*% t(C_xX)
  return(s)
}

#-----------------------------------------
#Predicted values
#-----------------------------------------

Vrv_pred = Vrv
Vlv_pred = Vlv
time_pred = tsteps

Pred_data = matrix(nrow = length(Vrv_pred), ncol = 4)

for (i in 1:length(Vrv_pred)){
  
  Pred_data[i,1] = pred_mean_SE(Vrv_pred[i], Vrv_obs, Vlv_pred[i], Vlv_obs, time_pred[i], tsteps_obs, Vspt_obs)
  
  Pred_data[i,2] = pred_cov_SE(Vrv_pred[i], Vrv_obs, Vlv_pred[i], Vlv_obs, time_pred[i], tsteps_obs)
  
  Pred_data[i,3] = Pred_data[i,1] + (1.997 * sqrt(Pred_data[i,2]))
  
  Pred_data[i,4] = Pred_data[i,1] - (1.997 * sqrt(Pred_data[i,2]))
}

SE_plot = ggplot() + 
  geom_line(aes(x = time_pred, y = Pred_data[,1]), col = "#007d69", lwd = 0.8) + 
  geom_point(aes(x = tsteps_obs, y = Vspt_obs), col = "black") + 
  geom_line(aes(x = time_pred, y = Vspt), col = "black", alpha = 0.5) +
  geom_line(aes(x = time_pred, y = Pred_data[,3]), lty = 2, col = "#813534") + 
  geom_line(aes(x = time_pred, y = Pred_data[,4]), lty = 2, col = "#813534") + 
  labs(title = "SE kernel interpolation",
       x = "Time",
       y = "Vspt")

SE_plot

mean_vector = pred_mean_QP(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs, Vspt_obs)
cov_mat = pred_cov_QP(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs)

posterior_samples = TruncatedNormal::rtmvnorm(n = 5, mu = mean_vector, sigma = Cov_matrix,
                                              lb = rep(0, length(Mean_vector)), 
                                              ub = rep(Inf, length(Mean_vector)))
