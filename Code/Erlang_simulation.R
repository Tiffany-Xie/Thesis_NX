# Erlang Distribution
# Packages ####
library(viridis)
library(ggplot2)
library(deSolve)


# Exact Distribution ####
erlang <- function(x, n, γ) {
  (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

plot_erlang <- function(x, n, γ) {
  if (length(n) > 1 & length(γ) > 1) {
    return("Only one parameter variable everytime :)")
  }
  df <- expand.grid(Time = x, Shape = n, Rate = γ)
  df$P <- with(df, erlang(Time, Shape, Rate))
  par_n <- max(length(n), length(γ))
  rand_col <- viridis(par_n)
  if (length(n) > length(γ)) (
    ggplot(df, aes(x=Time, y=P, color=factor(Shape))) + geom_line() +
      scale_color_manual(values = rand_col) +
      geom_vline(xintercept = 1/γ, color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = paste("Probability density functions for Erlang distributions with mean (1/γ) =", 1/γ),
           x = "Time (days)",
           y = "Probability Density",
           color = "Shape Parameter (n)")
  )
  else if (length(γ) > length(n)) {
    ggplot(df, aes(x=Time, y=P, color=factor(Rate))) + geom_line() +
      scale_color_manual(values = rand_col) +
      labs(title = paste("Probability density functions for Erlang distributions with shape (n) =", n),
           x = "Time (days)",
           y = "Probability Density",
           color = "Scale Parameter (γ)")
  }
  else {
    ggplot(df, aes(x=Time, y=P)) + geom_line() +
      labs(title = paste("Probability density functions for Erlang distributions with shape (n) =", n, "and mean (1/γ) =", 1/γ),
           x = "Time (days)",
           y = "Probability Density")
  }
}

# Parameters
x <- seq(0,30, by=0.1)
n <- c(1,2,3,4,5,10,20,50,100)
γ <- 0.077
#n <- 10
#γ <- seq(0.1, 0.7, by=0.1)
#n <- 2
#γ <- 0.1
plot_erlang(x, n, γ)


# Exact vs. Expact ####
model <- function(t, states, params) {
  with(as.list(c(states, params)), {
    dI1 <- -(10*γ + μ)*I1
    dI2 <- 10*γ*I1 - (10*γ + μ)*I2
    dI3 <- 10*γ*I2 - (10*γ + μ)*I3
    dI4 <- 10*γ*I3 - (10*γ + μ)*I4
    dI5 <- 10*γ*I4 - (10*γ + μ)*I5
    dI6 <- 10*γ*I5 - (10*γ + μ)*I6
    dI7 <- 10*γ*I6 - (10*γ + μ)*I7
    dI8 <- 10*γ*I7 - (10*γ + μ)*I8
    dI9 <- 10*γ*I8 - (10*γ + μ)*I9
    dI10 <- 10*γ*I9 - (10*γ + μ)*I10
    dR <- 10*γ*I10 - μ*R
    dR_delta <- dR - R_delta
    R_delta <- dR
    list(c(dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dR, dR_delta)) # dR_delta 
  })
}

compare <- function(time, n, γ, μ, states, model) {
  df <- data.frame(Time = time)
  df$P <- erlang(df$Time, n, γ)
  soln <- ode(y = states,
              times = time, 
              func = model, 
              parms = c(γ, μ))
  df$P_est <- soln[,"R_delta"]
  ggplot(df, aes(x=Time)) +
    geom_line(aes(y = P, col = "Exact P"), linewidth = 1) +
    geom_line(aes(y = P_est, col = "Estimated P"), linewidth = 1, linetype = "twodash") +
    labs(title = paste("Probability density functions with 1/γ =", 1/γ, "and n =", n), 
         x = "Stage duration (days)",
         y = "Probability Density")
}

states <- c(I1 = 1, I2 = 0, I3 = 0, 
            I4 = 0, I5 = 0, I6 = 0,
            I7 = 0, I8 = 0, I9 = 0,
            I10 = 0, R = 0, R_delta = 0)
time <- seq(0,30,by=0.1)
n <- 10
γ <- 0.1
μ <- 0
compare(time, n, γ, μ, states, model)





