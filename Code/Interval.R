library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)

######################################################################
## gamma -> gamma
######################################################################

set.seed(223)

numObs <- 1000
kappa <- 0.3
D <- 7

## k = 1/shape
## D = shape*scale = scale/kappa

v <- round(rgamma(numObs, shape=1/kappa, scale=D*kappa))

print(mean(v))
print(sd(v)^2/mean(v)^2) # nice!

## start/end times
## start is a uniform
## end is start + v [above]
## round them both before subtracting

######################################################################

gammaDelay <- function(n, mean, kappa){
  set.seed(223)
  v <- round(rgamma(n, shape=1/kappa, scale=mean*kappa))
  print(c(mean = mean, kappa = kappa,
        Dmean = mean(v), Dkappa = sd(v)^2/mean(v)^2))
  return(v)
}

gg.nll <- function(logmean, logkappa, interval) {
  -sum(dgamma(interval, shape = 1/exp(logkappa), scale = exp(logmean+logkappa), log=TRUE))
}

mean = 7
kappa = 0.3
n = 1000
g_interval = gammaDelay(n,mean,kappa) # v
startPar = list(logmean = 2, logkappa = -1)
fit <- mle2(gg.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

time = seq(min(g_interval), max(g_interval), length.out=n)
fitmean = exp(coef(fit)[["logmean"]])
fitkappa = exp(coef(fit)[["logkappa"]])
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa)*n,
                 fit_gamma = dgamma(time, shape=1/fitkappa, scale=fitmean*fitkappa)*n)
ggplot(df) + 
  geom_histogram(aes(x=interval)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Gamma after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> gamma")

######################################################################
## gamma -> Erlang
######################################################################








