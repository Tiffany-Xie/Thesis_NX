# Test Pesudo Erlang PDF
library(pseudoErlang)
library(ggplot2); theme_set(theme_minimal(base_size = 13))

erlang <- function(x, n, γ) {
  (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

β <- 0
D <- 10
n <- 4
kappa <- 0.5
μ <- 0.00
N <- 1
I0 <- 1
ts <- 1
T <- 100

perlang <- function(t, n, D, kappa) {
  r <- kappa2r(kappa, n)
  a <- (1-1/r^n)/(D*(1-1/r))
  ans <- 0
  
  for (i in 1:n) {
    innerp <- 1
    
    for (j in 1:n) {
      if (j != i) {
        innerp <- innerp * r^(j-1)/(r^(j-1) - r^(i-1))
      }
    }
    ans <- ans + innerp * a * r^(i-1) * exp(-t*a*r^(i-1))
  }
  return(ans)
}

tPEpdf <- function(D, kappa, n, ts, T) {
  sigr <- sinnerFlow(β=0, D=D, kappa=kappa, n=n, μ=0, N=1, I0=1, ts=ts, T=T)
  time <- timeSeq(ts, T)
  PEpdf <- perlang(time, n, D, kappa)
  df <- data.frame(Time = time,
                   ODE = diff(sigr[,"R"]),
                   PEpdf = PEpdf)
  print(ggplot(df, aes(x=Time)) +
          geom_line(aes(y=ODE, color = "ODE"), linewidth=1) +
          geom_line(aes(y=PEpdf, color = "PE PDF"), linewidth=2, linetype="dashed") +
          labs(y = "Probability Density")
        )
  return(df)
}

tPEpdf(D, kappa, n, ts, T)
tPEpdf(D, kappa, n=5, ts, T)
tPEpdf(D, kappa=0.2, n=7, ts, T)
#ggplot(df, aes(x=Time, y=PEpdf)) + geom_line()
df = tPEpdf(D, kappa=0.2, n=7, ts, T)
