# CVS system of ordinary differential equations (ODEs)
# load all the necessary packages
library(deSolve)
library(ggplot2)
library(tidyverse)
library(patchwork)

theme_set(theme_bw())

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
pred_mean_QP <- function(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs, y, nugget = tau_sq_MAP){
  C_pred_obs <- Cov_QP(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs)
  C_obs_obs <- Cov_QP(x1_obs, x1_obs, x2_obs, x2_obs, x3_obs, x3_obs)
  m <-C_pred_obs %*% solve(C_obs_obs + (nugget*diag(ncol(C_obs_obs)))) %*% y
  return(m)
}

#Posterior covariance matrix using Cholesky Decomposition
pred_cov_QP <- function(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs, 
                        nugget = tau_sq_MAP) {
  
  C_pred_obs <- Cov_QP(x1_pred, x1_obs, x2_pred, x2_obs, x3_pred, x3_obs)
  C_obs_obs  <- Cov_QP(x1_obs, x1_obs, x2_obs, x2_obs, x3_obs, x3_obs)
  C_pred_pred <- Cov_QP(x1_pred, x1_pred, x2_pred, x2_pred, x3_pred, x3_pred)
  
  # Add nugget for numerical stability
  C_obs_obs_nug <- C_obs_obs + nugget * diag(nrow(C_obs_obs))
  
  # Cholesky decomposition for inversion
  L <- chol(C_obs_obs_nug)
  tmp <- forwardsolve(t(L), t(C_pred_obs))  # Solving Lᵗ x = Kᵗ
  s <- C_pred_pred - t(tmp) %*% tmp
  
  return(s)
}

# Smith cardiovascular model:
CVS <- function(time, u, params) {
  
  # unpack the state variables
  QmtQP <- u[1]; QavQP <- u[2]; QtcQP <- u[3]; QpvQP <- u[4]
  VlvQP <- u[5]; VaoQP <- u[6]; VvcQP <- u[7]; VrvQP <- u[8]
  VpaQP <- u[9]; VpuQP <- u[10]
  
  # Unpack the parameters
  with(as.list(params), {
    
    # Activation function
    e <- exp(-80 * (time %% 0.75 - 0.375)^2)
    
    # Pericardial pressure calculation
    Vpcd <- VlvQP + VrvQP
    Ppcd <- P0pcd * (exp(lambdapcd * (Vpcd - V0pcd)) - 1)
    Pperi <- Ppcd + Pth
    
    # Solve for septum volume using GP
    mean_conditional = pred_mean_QP(VrvQP,Vrv_obs,VlvQP,Vlv_obs,time,tsteps_obs,Vspt_obs)
    cov_conditional = pred_cov_QP(VrvQP,Vrv_obs,VlvQP,Vlv_obs,time,tsteps_obs)
    
    VsptQP = mean_conditional + (3 * sqrt(cov_conditional))
    
    # Calculating filling volumes
    VlvQPf <- VlvQP - VsptQP
    VrvQPf <- VrvQP + VsptQP
    
    # Calculating pressures
    Plvf <- e * Elvf * (VlvQPf - Vdlvf) + (1 - e) * P0lvf * (exp(lambdalvf * (VlvQPf - V0lvf)) - 1)
    Prvf <- e * Ervf * (VrvQPf - Vdrvf) + (1 - e) * P0rvf * (exp(lambdarvf * (VrvQPf - V0rvf)) - 1)
    
    Plv <- Plvf + Pperi
    Prv <- Prvf + Pperi
    
    Pao <- Eao * (VaoQP - Vdao)
    Pvc <- Evc * (VvcQP - Vdvc)
    Ppa <- Epa * (VpaQP - Vdpa) + Pth
    Ppu <- Epu * (VpuQP - Vdpu) + Pth
    
    Qsys <- (Pao - Pvc) / Rsys
    Qpul <- (Ppa - Ppu) / Rpul
    
    # Differential equations with conditions
    du <- numeric(10)
    du[1] <- ifelse(Ppu - Plv > 0 || QmtQP > 0, (Ppu - Plv - (QmtQP * Rmt)) / Lmt, 0)
    du[2] <- ifelse(Plv - Pao > 0 || QavQP > 0, (Plv - Pao - (QavQP * Rav)) / Lav, 0)
    du[3] <- ifelse(Pvc - Prv > 0 || QtcQP > 0, (Pvc - Prv - (QtcQP * Rtc)) / Ltc, 0)
    du[4] <- ifelse(Prv - Ppa > 0 || QpvQP > 0, (Prv - Ppa - (QpvQP * Rpv)) / Lpv, 0)
    
    # Flow adjustments
    QmtQP <- max(QmtQP, 0)
    QavQP <- max(QavQP, 0)
    QtcQP <- max(QtcQP, 0)
    QpvQP <- max(QpvQP, 0)
    
    # Volume differential equations
    du[5] <- QmtQP - QavQP
    du[6] <- QavQP - Qsys
    du[7] <- Qsys - QtcQP
    du[8] <- QtcQP - QpvQP
    du[9] <- QpvQP - Qpul
    du[10] <- Qpul - QmtQP
    
    return(list(du))
  })
}


# Define time range
tspan <- c(0.0, 10)
tsteps <- seq(tspan[1], tspan[2], length.out = 1000)

# Initial conditions
u0 <- c(QmtQP = 245.5813, QavQP = 0, QtcQP = 190.0661, QpvQP = 0, VlvQP = 94.6812, VaoQP = 133.3381, 
        VvcQP = 329.7803, VrvQP = 90.7302, VpaQP = 43.0123, VpuQP = 808.4579)

# Parameters
p_ <- list(Elvf = 2.8798, Eao = 0.6913, Evc = 0.0059, Ervf = 0.585, Epa = 0.369, 
           Epu = 0.0073, Rmt = 0.0158, Rav = 0.018, Rsys = 1.0889, Rtc = 0.0237, 
           Rpv = 0.0055, Rpul = 0.1552, Lmt = 7.6968e-5, Lav = 1.2189e-4, Ltc = 8.0093e-5, 
           Lpv = 1.4868e-4, Vdlvf = 0, Vdao = 0, Vdvc = 0, Vdrvf = 0, Vdpa = 0, Vdpu = 0, 
           P0lvf = 0.1203, P0rvf = 0.2157, lambdalvf = 0.033, lambdarvf = 0.023, 
           Espt = 48.754, V0lvf = 0, V0rvf = 0, P0spt = 1.1101, P0pcd = 0.5003, 
           V0spt = 2, V0pcd = 200, lambdaspt = 0.435, lambdapcd = 0.03, Vdspt = 2, 
           Pth = -4)


# Solve the system
system.time(sol <- ode(y = u0, times = tsteps, func = CVS, parms = p_, method = "ode45"))
DataQPp3 <- as.data.frame(sol)

# Vspt and calculation
QmtQPp3 <- DataQPp3$QmtQP
QavQPp3 <- DataQPp3$QavQP
QtcQPp3 <- DataQPp3$QtcQP
QpvQPp3 <- DataQPp3$QpvQP
VlvQPp3 <- DataQPp3$VlvQP
VaoQPp3 <- DataQPp3$VaoQP
VvcQPp3 <- DataQPp3$VvcQP
VrvQPp3 <- DataQPp3$VrvQP
VpaQPp3 <- DataQPp3$VpaQP
VpuQPp3 <- DataQPp3$VpuQP

VsptQPp3 <- numeric(length(tsteps))

for (i in seq_along(tsteps)) {
  mean_conditional = pred_mean_QP(VrvQPp3[i],Vrv_obs,VlvQPp3[i],Vlv_obs,tsteps[i],tsteps_obs,Vspt_obs)
  cov_conditional = pred_cov_QP(VrvQPp3[i],Vrv_obs,VlvQPp3[i],Vlv_obs,tsteps[i],tsteps_obs)
  
  VsptQPp3[i] = mean_conditional + (3 * sqrt(cov_conditional))
}

df_Vspt <- tibble(time = tsteps, 
                  VlvQP = VlvQP, 
                  VrvQP = VrvQP, 
                  VaoQP = VaoQP, 
                  VsptQP = VsptQP)

ggplot() + 
  geom_line(aes(x = tsteps, y = Vspt, color = "Newton Raphson", linetype = "Newton Raphson")) + 
  geom_line(aes(x = tsteps, y = VsptQPp3, color = "QP Kernel", linetype = "QP Kernel"), lwd = 1.5) + 
  geom_line(aes(x = tsteps, y = VsptQPp3, color = "SE Kernel", linetype = "SE Kernel"), lwd = 1) + 
  labs(title = "Temporal evolution of Vspt", x = "Time (s)", y = "Vspt (ml)", color = "Method", linetype = "Method") +
  scale_color_manual(values = c("Newton Raphson" = "black", "QP Kernel" = "red", "SE Kernel" = "blue")) +
  scale_linetype_manual(values = c("Newton Raphson" = "solid", "QP Kernel" = "dotted", "SE Kernel" = "longdash")) +
  theme_minimal()

VsptQPp3

VlvQPp3
