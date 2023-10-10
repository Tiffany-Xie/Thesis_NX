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
#n <- 3
#γ <- 0.1
plot_erlang(x, n, γ)

#newmodel <- function(t, states, params) {
#  with(as.list(c(params)), {
#    boxes <- length(states) - 2
#    I <- states[1:boxes]
#    R <- states[boxes+1]
#    C <- states[boxes+2]
    
#    dI <- -(boxes*x+μ)*I
#    dI[2:boxes] <- dI[2:boxes]
#  })
#}

# Exact vs. Expact ####
model <- function(t, states, params) {
  with(as.list(c(states, params)), {
    dI <- numeric(n)
    #states <- c(states[1], numeric(n-1), states[2])
    dI[1] <- -(n*γ + μ)*states[1]
    if (n>1) {
      for (i in 2:n) {
        dI[i] <- n*γ*states[i-1] - (n*γ + μ)*states[i]
      }
    }
    dR <- n*γ*states[n] - μ*R
    list(c(dI, dR))
  })
}

compare <- function(time, n, γ, μ, states, model) {
  params <- c(γ = γ, μ = μ, n = n)
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
    labs(title = paste("Probability density functions with 1/γ =", 1/γ, "and n =", n), 
         x = "Stage duration (days)",
         y = "Probability Density")
}


states <- c(I1 = 1, I2 = 0, I3 = 0, 
            I4 = 0, I5 = 0, I6 = 0,
            I7 = 0, R = 0)
time <- seq(0,30,by=0.01)
n <- 7
γ <- 0.1
μ <- 0
compare(time, n, γ, μ, states, model)


# Gamma vs. Erlang (check)
#library(tidyr)

times <- seq(0,30,by=0.01)
n <- 3
γ <- 0.1
df <- data.frame(Time = times,
                 Gamma = dgamma(times, n, n*γ),
                 Erlang = erlang(times, n, γ))
#df <- tidyr::gather(df, key = "Type", value = "PDF", -Time)
#ggplot(df, aes(x = Time, y = PDF, color = Type)) +
#  geom_line() +
#  labs(title = paste("Gamma vs. Erlang, with shape(n) =", n, "and rate(γ) =", γ), 
#       x = "Time (days)",
#       y = "Probability Density") +
#  theme_minimal()

ggplot(df, aes(x=Time)) +
  geom_line(aes(y = Gamma, col = "dgamma"), linewidth = 1) +
  geom_line(aes(y = Erlang, col = "formula"), linewidth = 2, linetype = "twodash") +
  labs(title = paste("Gamma vs. Erlang, with shape(n) =", n, "and rate(γ) =", γ), 
       x = "Stage duration (days)",
       y = "Probability Density")

print(
	ggplot(df) + aes(x=Time)
	+ geom_line(aes(y = Gamma, col = "dgamma"), linewidth = 1)
	+ geom_line(aes(y = Erlang, col = "formula")
		, linewidth = 2, linetype = "twodash"
	)
	+ labs(
		title = paste("Gamma vs. Erlang, with shape(n) =", n, "and rate(γ) =", γ) 
		, x = "Stage duration (days)"
		, y = "Probability Density"
	)
)


