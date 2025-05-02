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
mcmc_hist(param_samples$draws(variables = c("sigma", "length_scale1", "length_scale2","length_scale3")))

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
  
  C_pred_obs <- Cov_Exp_sq(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs)
  C_obs_obs  <- Cov_Exp_sq(x1_obs, x1_obs, x2_obs, x2_obs, x3_obs, x3_obs)
  C_pred_pred <- Cov_Exp_sq(x1_pred, x1_pred, x2_pred, x2_pred, x3_pred, x3_pred)
  
  # Add nugget for numerical stability
  C_obs_obs_nug <- C_obs_obs + nugget * diag(nrow(C_obs_obs))
  
  # Cholesky decomposition for inversion
  L <- chol(C_obs_obs_nug)
  tmp <- forwardsolve(t(L), t(C_pred_obs))
  s <- C_pred_pred - t(tmp) %*% tmp
  
  return(s)
}

#-----------------------------------------
#Predicted values
#-----------------------------------------

Vrv_pred = Vrv
Vlv_pred = Vlv
time_pred = tsteps

Pred_data = matrix(nrow = length(Vrv_pred), ncol = 4)

Pred_data[,1] = pred_mean_SE(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs, Vspt_obs)

Pred_data[,2] = diag(pred_cov_SE(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs))

Pred_data[,3] = Pred_data[,1] + (3 * sqrt(Pred_data[,2]))

Pred_data[,4] = Pred_data[,1] - (3 * sqrt(Pred_data[,2]))

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

SEMeanMSE = mean((Pred_data[,1] - Vspt)^2)

#Compute mean vector and cov matrix for MVN
Cov_matrix = pred_cov_SE(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs) + 1e-7 * diag(nrow(Pred_df))
Mean_vector = pred_mean_SE(Vrv_pred, Vrv_obs, Vlv_pred, Vlv_obs, time_pred, tsteps_obs, Vspt_obs)

#sample

Vspt_samplesSE = list()
for (i in 1:10){
  Vspt_samplesSE[[i]] = rmvnorm(n = 1, mu = Mean_vector, Sigma = Cov_matrix)
}

Vspt_samples_mat = do.call(rbind, Vspt_samplesSE)

SEsampleMSE = numeric(10)

for (i in 1:10){
  
  SEsampleMSE[i] = mean((Vspt - Vspt_samples_mat[i,])^2)
  
}

SEsampleMSE

mean(SEsampleMSE)

plot(tsteps, Vspt_samples_mat[1,], type = "l")
for (i in 2:10){
  lines(tsteps, Vspt_samples_mat[i,], type = "l")
}

