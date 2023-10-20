# Packages ####
library(viridis)
library(ggplot2)
library(deSolve)
library(gridExtra)

library(shellpipes)
startGraphics(width=11)

# Models ####
erlang <- function(x, n, γ) {
  (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

SInR_geom <- function(t, states, params) {
  with(as.list(c(params)), {
    I <- states[1:n]
    R <- states[[n+1]]
    
    Iprev <- c(0, I[1:(n-1)])
    dI <- a*r^(c(0,0:(n-2)))*Iprev - (a*r^(0:(n-1))+μ)*I
    dR <- a*r^(n-1)*I[[n]] - μ*R
    return(list(c(dI, dR)))
  })
}

# Compare ####
EP_compare <- function(time, γ, μ, r, nE, nPE, model) {
  df_E <- data.frame(Time = time)
  df_E$P_E <- erlang(time, nE, γ) 
  
  df_PE <- expand.grid(Time = time, nPE= nPE, nE = nE, γ = γ, r = r)
  df_PE$a <- with(df_PE, (1/r^nPE-1)*γ/(1/r-1)) #sqrt(nE*γ^2*((1/r^2)^nPE - 1)/(1/r^2-1))
  states <- c(1, numeric(nPE))
  names(states) <- c(paste ("I", 1:nPE, sep=""), "R")
  P_PE <- c()
  for (i in r) {
    a <- (1/i^nPE-1)*γ/(1/i-1)
    params <- c(γ = γ, μ = μ, n = nPE, a = a, r = i)
    soln <- ode(y = states,
                times = time, 
                func = model, 
                parms = params)
    P_PE <- c(P_PE, c(diff(soln[,"R"])/diff(time), NA))
  }

  df_PE$P_PE <- P_PE 
  
  p1 <- ggplot(df_PE, aes(x=Time, y=P_PE, color=factor(r))) +
    geom_line() +
    scale_color_manual(values = viridis(length(r))) +
    geom_vline(xintercept = 1/γ, color = "red", linetype = "dashed", linewidth = 1) +
    geom_line(data=df_E, aes(x = Time, y = P_E), color = "black", linewidth = 1.5) +
    labs(title = paste0("Compare different r: D=", 1/γ, ", nE=", nE, ", nPE=", nPE),
         x = "Stage duration (days)",
         y = "Probability Density")
  
  df_r <- data.frame(r = r)
  df_r$a <- with(df_r, (1/r^nPE-1)*γ/(1/r-1))
  df_r$var <- with(df_r, 1/a^2*(1/r^(2*nPE)-1)/(1/r^2-1))
  df_r$mean <- with(df_r, 1/a*(1/r^nPE-1)/(1/r-1))
  df_r$K <- with(df_r, var/mean^2)
  p2 <- ggplot(df_r, aes(x=r, y=K)) + geom_line() +
    geom_hline(yintercept = 1/(nE), color = "red", linetype = "dashed", linewidth = 1) +
    labs(title="Coefficient of Variation^2 (K)")
  
  plots <- list()
  plots[[1]] <- p1
  plots[[2]] <- p2
  grid.arrange(grobs = plots, ncol = 2)
  
}

# Simulation ####
time = seq(0,30,by=0.01) 
γ <- 0.1 
μ <- 0
nPE <- 12

# nE == 2
r <- seq(2.5, 3.5, by=0.1)
nE <- 2
EP_compare(time, γ, μ, r, nE, nPE, SInR_geom)

# nE == 4
r <- seq(1.5, 2.5, by=0.1)
nE <- 4
EP_compare(time, γ, μ, r, nE, nPE, SInR_geom)

# nE == 8
r <- seq(1.2, 1.30, by=0.01)
nE <- 8
EP_compare(time, γ, μ, r, nE, nPE, SInR_geom)








