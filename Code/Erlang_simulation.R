# Packages ####
library(viridis)
library(ggplot2)
library(deSolve)

# Erlang Distribution ####
erlang <- function(x, n, γ) {
  (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

# Models & Compare ####
SInR <- function(t, states, params) {
  with(as.list(c(params)), {
    I <- states[1:n]
    R <- states[[n+1]]
    
    Iprev <- c(0, I[1:(n-1)])
    dI <- n*γ*Iprev - (n*γ+μ)*I
    dR <- n*γ*I[[n]] - μ*R
    return(list(c(dI, dR)))
  })
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

model_geom <- function(t, states, params) {
  with(as.list(c(states, params)), {
    dI <- numeric(n)
    dI[1] <- -(a + μ)*states[1]
    if (n>1) {
      for (i in 2:n) {
        dI[i] <- a*r^(i-2)*states[i-1] - (a*r^(i-1) + μ)*states[i]
      }
    }
    dR <- a*r^(n-1)*states[n] - μ*R
    list(c(dI, dR))
  })
}

compare <- function(time, n, γ, μ, a, r, states, model) {
  params <- c(γ = γ, μ = μ, n = n, a = a, r = r)
  df <- data.frame(Time = time)
  df$P <- erlang(time, n, γ)
  soln <- ode(y = states,
              times = time, 
              func = model, 
              parms = params)
  
  df$P_est <- c(diff(soln[,"R"])/diff(time), NA)
  ggplot(df, aes(x=Time)) +
    geom_line(aes(y = P, col = "Exact P"), linewidth = 1) +
    geom_line(aes(y = P_est, col = "ODE"), linewidth = 2, linetype = "twodash") +
    labs(title = paste("Probability density functions with D =", 1/γ, "and stage =", n), 
         x = "Stage duration (days)",
         y = "Probability Density")
}


# Geometric γ ####
time <- seq(0,30,by=0.01)
n <- 4
γ <- 0.1
μ <- 0
states <- c(1, numeric(n))
names(states) <- c(
  paste ("I", 1:n, sep="")
  , "R"
)

# Geometric: r = 2, a = 3/16, γ_i ascending = TRUE, D = 10, n = 4
r <- 2
a <- 3/16 # 127/640
compare(time, n, γ, μ, a, r, states, SInR_geom)

# Geometric: different r values, to maintain D, a = (1/r^n-1)*γ/(1/r-1), D = 10, n =  4
r <- seq(1.2,3,by=0.2)
df <- expand.grid(Time = time, n = n, γ = γ, r = r)
df$a <- with(df, (1/r^n-1)*γ/(1/r-1)) 
df$P <- with(df, erlang(Time, n, γ))
ODE <- c()
for (i in r) {
  soln <- ode(y = states,
              times = time, 
              func = SInR_geom, 
              parms = c(γ = γ, μ = μ, n = n, a = (1/i^n-1)*γ/(1/i-1), r = i))
  ODE <- c(ODE, c(diff(soln[,"R"])/diff(time), NA))
}
df$ODE <- ODE

ggplot(df, aes(x = Time, y = ODE, color = factor(r))) +
  geom_line() +
  scale_color_manual(values = viridis(length(r))) +
  geom_vline(xintercept = 1/γ, color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(aes(x = Time, y = P), color = "black", linewidth = 1.5) +
  labs(title = paste0("Compare different r: D=", 1/γ[1], ", n=", n[1]),
       x = "Stage duration (days)",
       y = "Probability Density")

# Geometric: same r, D; different n
library(gridExtra)
N <- seq(2,7)
r <- 2
plots <- list()
i = 1
for (n in N) {
  a = (1/r^n-1)*γ/(1/r-1)
  states <- c(1, numeric(n))
  names(states) <- c(
    paste ("I", 1:n, sep="")
    , "R"
  )
  p <- compare(time, n, γ, μ, a, r, states, SInR_geom)
  plots[[i]] <- p
  i <- i+1
}
grid.arrange(grobs = plots, ncol = 2)


