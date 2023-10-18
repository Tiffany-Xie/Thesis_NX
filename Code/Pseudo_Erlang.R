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
  df_E <- data.frame(Time = time)
  df_E$P_E <- erlang(time, nE, γ)
  
  df_SE <- expand.grid(Time = time, nSE= nSE, nE = nE, γ = γ, r = r)
  df_SE$a <- with(df_SE, sqrt(nE*γ^2*((1/r^2)^nSE - 1)/(1/r^2-1)))
  states <- c(1, numeric(nSE))
  names(states) <- c(paste ("I", 1:nSE, sep=""), "R")
  P_SE <- c()
  for (i in r) {
    a <- sqrt(nE*γ^2*((1/i^2)^nSE - 1)/(1/i^2-1))
    params <- c(γ = γ, μ = μ, n = nSE, a = a, r = i)
    soln <- ode(y = states,
                times = time, 
                func = model, 
                parms = params)
    P_SE <- c(P_SE, c(diff(soln[,"R"])/diff(time), NA))
  }

  df_SE$P_SE <- P_SE
  
  ggplot(df_SE, aes(x=Time, y=P_SE, color=factor(r))) +
    geom_line() +
    scale_color_manual(values = viridis(length(r))) +
    geom_vline(xintercept = 1/γ, color = "red", linetype = "dashed", linewidth = 1) +
    geom_line(data=df_E, aes(x = Time, y = P_E), color = "black", linewidth = 1.5) +
    labs(title = paste0("Compare different r: D=", 1/γ, ", nE=", nE, ", nSE=", nSE),
         x = "Stage duration (days)",
         y = "Probability Density")
}

# Simulation ####
time = seq(0,30,by=0.01) 
γ <- 0.1 
μ <- 0
#r <- 2
r <- c(1.5, 2, 2.5, 3)
nE <- 4
nSE <- 12

ES_compare(time, γ, μ, r, nE, nSE, SInR_geom)




#### Try a,r variable cal
# Example nonlinear system
equations <- function(x) {
  y1 <- x[1]^2 + x[2]^2 - 25
  y2 <- x[1]*x[2] - 9
  return(c(y1, y2))
}

# Convert the equations to a formula
formula <- as.formula(c(y1, y2) ~ x[1] + x[2])

# Initial guess
initial_guess <- c(x1 = 0, x2 = 0)

quit()

# Solve the system of nonlinear equations
solution <- nls(formula, start = initial_guess)

# Print the result
print(summary(solution))
