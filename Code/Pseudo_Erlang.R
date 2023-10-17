# Packages ####
library(viridis)
library(ggplot2)
library(deSolve)

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
ES_compare <- function(time, γ, μ, r, nE, nSE, model) {
  a <- sqrt(nE*γ^2*((1/r^2)^nSE - 1)/(1/r^2-1))
  params <- c(γ = γ, μ = μ, n = nSE, a = a, r = r)
  states <- c(1, numeric(nSE))
  names(states) <- c(
    paste ("I", 1:nSE, sep="")
    , "R"
  )
  
  df <- data.frame(Time = time)
  df$P_E <- erlang(time, nE, γ)
  
  soln <- ode(y = states,
              times = time, 
              func = model, 
              parms = params)
  df$P_SE <- c(diff(soln[,"R"])/diff(time), NA)
  
  ggplot(df, aes(x=Time)) +
    geom_line(aes(y = P_E, col = "Erlang"), linewidth = 1) +
    geom_line(aes(y = P_SE, col = "Pseudo Erlang"), linewidth = 2, linetype = "twodash") +
    labs(title = paste0("Probability density functions with D=", 1/γ, ", Erlang stage=", nE, ", Pseudo stage=", nSE), 
         x = "Stage duration (days)",
         y = "Probability Density")
}

# Simulation ####
time = seq(0,30,by=0.01) 
γ <- 0.1 
μ <- 0
r <- 2
nE <- 4
nSE <- 12

ES_compare(time, γ, μ, r, nE, nSE, SInR_geom)





