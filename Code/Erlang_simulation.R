# Erlang Distribution
library(viridis)
library(ggplot2)

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

