# CVS system of ordinary differential equations (ODEs)
# load all the necessary packages
library(deSolve)
library(ggplot2)
library(tidyverse)
library(patchwork)

theme_set(theme_bw())

# Smith cardiovascular model:
CVS <- function(time, u, params) {
  
  # unpack the state variables
  Qmt <- u[1]; Qav <- u[2]; Qtc <- u[3]; Qpv <- u[4]
  Vlv <- u[5]; Vao <- u[6]; Vvc <- u[7]; Vrv <- u[8]
  Vpa <- u[9]; Vpu <- u[10]
  
  # Unpack the parameters
  with(as.list(params), {
    
    # Activation function
    e <- exp(-80 * (time %% 0.75 - 0.375)^2)
    
    # Pericardial pressure calculation
    Vpcd <- Vlv + Vrv
    Ppcd <- P0pcd * (exp(lambdapcd * (Vpcd - V0pcd)) - 1)
    Pperi <- Ppcd + Pth
    
    # Solve for septum volume using Newton-Raphson method
    Vspt <- 0
    dy <- 0.01
    while (abs(dy) > 1e-5 * abs(Vspt)) {
      f <- e * Espt * (Vspt - Vdspt) + (1 - e) * P0spt * (exp(lambdaspt * (Vspt - V0spt)) - 1) -
        e * Elvf * (Vlv - Vspt - Vdlvf) - (1 - e) * P0lvf * (exp(lambdalvf * (Vlv - Vspt - V0lvf)) - 1) +
        e * Ervf * (Vrv + Vspt - Vdrvf) + (1 - e) * P0rvf * (exp(lambdarvf * (Vrv + Vspt - V0rvf)) - 1)
      
      df <- e * Espt + lambdaspt * (1 - e) * P0spt * exp(lambdaspt * (Vspt - V0spt)) +
        e * Elvf + lambdalvf * (1 - e) * P0lvf * exp(lambdalvf * (Vlv - Vspt - V0lvf)) +
        e * Ervf + lambdarvf * (1 - e) * P0rvf * exp(lambdarvf * (Vrv + Vspt - V0rvf))
      
      dy <- f / df
      Vspt <- Vspt - dy
    }
    
    # Calculating filling volumes
    Vlvf <- Vlv - Vspt
    Vrvf <- Vrv + Vspt
    
    # Calculating pressures
    Plvf <- e * Elvf * (Vlvf - Vdlvf) + (1 - e) * P0lvf * (exp(lambdalvf * (Vlvf - V0lvf)) - 1)
    Prvf <- e * Ervf * (Vrvf - Vdrvf) + (1 - e) * P0rvf * (exp(lambdarvf * (Vrvf - V0rvf)) - 1)
    
    Plv <- Plvf + Pperi
    Prv <- Prvf + Pperi
    
    Pao <- Eao * (Vao - Vdao)
    Pvc <- Evc * (Vvc - Vdvc)
    Ppa <- Epa * (Vpa - Vdpa) + Pth
    Ppu <- Epu * (Vpu - Vdpu) + Pth
    
    Qsys <- (Pao - Pvc) / Rsys
    Qpul <- (Ppa - Ppu) / Rpul
    
    # Differential equations with conditions
    du <- numeric(10)
    du[1] <- ifelse(Ppu - Plv > 0 || Qmt > 0, (Ppu - Plv - (Qmt * Rmt)) / Lmt, 0)
    du[2] <- ifelse(Plv - Pao > 0 || Qav > 0, (Plv - Pao - (Qav * Rav)) / Lav, 0)
    du[3] <- ifelse(Pvc - Prv > 0 || Qtc > 0, (Pvc - Prv - (Qtc * Rtc)) / Ltc, 0)
    du[4] <- ifelse(Prv - Ppa > 0 || Qpv > 0, (Prv - Ppa - (Qpv * Rpv)) / Lpv, 0)
    
    # Flow adjustments
    Qmt <- max(Qmt, 0)
    Qav <- max(Qav, 0)
    Qtc <- max(Qtc, 0)
    Qpv <- max(Qpv, 0)
    
    # Volume differential equations
    du[5] <- Qmt - Qav
    du[6] <- Qav - Qsys
    du[7] <- Qsys - Qtc
    du[8] <- Qtc - Qpv
    du[9] <- Qpv - Qpul
    du[10] <- Qpul - Qmt
    
    return(list(du))
  })
}


# Define time range
tspan <- c(0.0, 10)
tsteps <- seq(tspan[1], tspan[2], length.out = 1000)

# Initial conditions
u0 <- c(Qmt = 245.5813, Qav = 0, Qtc = 190.0661, Qpv = 0, Vlv = 94.6812, Vao = 133.3381, 
        Vvc = 329.7803, Vrv = 90.7302, Vpa = 43.0123, Vpu = 808.4579)

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
sol <- ode(y = u0, times = tsteps, func = CVS, parms = p_, method = "ode45")
original_data <- as.data.frame(sol)

# Adding noise to data
noise_magnitude <- 0.0
sd_vals <- apply(original_data, 2, sd)
noisy_data <- as.data.frame(sweep(original_data, 2, noise_magnitude * sd_vals, "+") + 
                              noise_magnitude * matrix(rnorm(nrow(original_data) * ncol(original_data)), 
                                                       nrow = nrow(original_data), ncol = ncol(original_data)))


noisy_data_tb <- as_tibble(noisy_data) %>%
  pivot_longer(cols = Qmt:Vpu, 
               names_to = "variable", 
               values_to = "vals")

p_noisy_data_tb <- ggplot(noisy_data_tb) +
  geom_line(aes(x = time, y = vals), color = "grey") +
  geom_point(aes(x = time, y = vals)) +
  facet_grid(vars(variable), scales = "free")

p_noisy_data_tb


# Vspt and calculation
Qmt <- noisy_data$Qmt
Qav <- noisy_data$Qav
Qtc <- noisy_data$Qtc
Qpv <- noisy_data$Qpv
Vlv <- noisy_data$Vlv
Vao <- noisy_data$Vao
Vvc <- noisy_data$Vvc
Vrv <- noisy_data$Vrv
Vpa <- noisy_data$Vpa
Vpu <- noisy_data$Vpu

e <- exp(-80 * (tsteps %% 0.75 - 0.375)^2)
Vpcd <- Vlv + Vrv
Ppcd <- p_$P0pcd * (exp(p_$lambdapcd * (Vpcd - p_$V0pcd)) - 1)
Pperi <- Ppcd + p_$Pth

Vspt <- rep(0, length(tsteps))
dx <- rep(0.01, length(tsteps))

for (i in seq_along(tsteps)) {
  while (abs(dx[i]) > 1e-5 * abs(Vspt[i])) {
    f <- e[i] * p_$Espt * (Vspt[i] - p_$Vdspt) + (1 - e[i]) * p_$P0spt * (exp(p_$lambdaspt * (Vspt[i] - p_$V0spt)) - 1) -
      e[i] * p_$Elvf * (Vlv[i] - Vspt[i] - p_$Vdlvf) - (1 - e[i]) * p_$P0lvf * (exp(p_$lambdalvf * (Vlv[i] - Vspt[i] - p_$V0lvf)) - 1) +
      e[i] * p_$Ervf * (Vrv[i] + Vspt[i] - p_$Vdrvf) + (1 - e[i]) * p_$P0rvf * (exp(p_$lambdarvf * (Vrv[i] + Vspt[i] - p_$V0rvf)) - 1)
    
    df <- e[i] * p_$Espt + p_$lambdaspt * (1 - e[i]) * p_$P0spt * exp(p_$lambdaspt * (Vspt[i] - p_$V0spt)) +
      e[i] * p_$Elvf + p_$lambdalvf * (1 - e[i]) * p_$P0lvf * exp(p_$lambdalvf * (Vlv[i] - Vspt[i] - p_$V0lvf)) +
      e[i] * p_$Ervf + p_$lambdarvf * (1 - e[i]) * p_$P0rvf * exp(p_$lambdarvf * (Vrv[i] + Vspt[i] - p_$V0rvf))
    
    dx[i] <- f / df
    Vspt[i] <- Vspt[i] - dx[i]
  }
}

df_Vspt <- tibble(time = tsteps, 
                  Vlv = Vlv, 
                  Vrv = Vrv, 
                  Vao = Vao, 
                  Vspt = Vspt)

p1 <- ggplot(df_Vspt) +
  geom_line(aes(x = time, y = Vspt), color = "black", lwd = 0.8) + 
  labs(title = "Temporal evolution of Vspt",
       x = "Time (s)",
       y = "Septum free wall volume, Vspt (ml)")

p2a <- ggplot(df_Vspt) +
  geom_point(aes(x = Vlv, y = Vspt)) +
  labs(title = "Vspt plotted against Vlv",
       x = "Volume of left ventricle, Vlv (ml)",
       y = "Septum free wall volume, Vspt (ml)")

p2b <- ggplot(df_Vspt) +
  geom_line(aes(x = time, y = Vlv), color = "black", lwd = 0.8) +
  labs(title = "Temporal evolution of Vlv",
       x = "Time (s)",
       y = "Volume of left ventricle, Vlv (ml)")

p3a <- ggplot(df_Vspt) +
  geom_point(aes(x = Vrv, y = Vspt)) + 
  labs(title = "Vspt plotted against Vrv",
       x = "Volume of right ventricle, Vrv (ml)",
       y = "Septum free wall volume, Vspt (ml)")

p3b <- ggplot(df_Vspt) +
  geom_line(aes(x = time, y = Vrv), color = "black", lwd = 0.8) + 
  labs(title = "Temporal evolution of Vrv",
       x = "Time (s)",
       y = "Volume of right ventricle, Vrv (ml)")

p4 <- ggplot(df_Vspt) +
  geom_point(aes(x = Vao, y = Vspt))


p3a

