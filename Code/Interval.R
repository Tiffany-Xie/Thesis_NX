library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)

######################################################################
## gamma -> gamma
######################################################################

set.seed(223)

numObs <- 5000
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
n = 2000
g_interval = gammaDelay(n,mean,kappa) # v
startPar = list(logmean = 2.5, logkappa = -0.5)
fit <- mle2(gg.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

time = seq(0, max(g_interval), length.out=n)
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
derlang <- function(x, shape, rate, log=FALSE) { # shape = n; rate = lambda
  density <- rate^shape * x^(shape-1) * exp(-rate*x) / factorial(shape-1)
  if (log) {
    density <- log(density)
  }
  return(density)
}

ge.nll <- function(logmean, logkappa, interval) {
  -sum(derlang(interval, shape = round(1/exp(logkappa)), rate = 1/exp(logmean)*round(1/exp(logkappa)), log=TRUE))
}

startPar = list(logmean = 2.5, logkappa = -2) # 2 -1
fit <- mle2(ge.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = exp(coef(fit)[["logmean"]])
fitkappa = exp(coef(fit)[["logkappa"]])
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa)*n,
                 fit_gamma = derlang(time, shape=round(1/fitkappa), rate=1/fitmean*round(1/fitkappa))*n)
ggplot(df) + 
  geom_histogram(aes(x=interval)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Erlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> Erlang")

#quit()
######################################################################
## gamma -> Erlang (2)
######################################################################
derlang2 <- function(x, D, shape=4, log=FALSE) { # shape = n; rate = lambda
  rate <- shape * 1/D
  density <- rate^shape * x^(shape-1) * exp(-rate*x) / factorial(shape-1)
  if (log) {
    density <- log(density)
  }
  return(density)
}

ge.nll2 <- function(logmean, interval) {
  -sum(derlang2(interval, D = exp(logmean), log=TRUE))
}

startPar = list(logmean = 1) # 2 -1
fit <- mle2(ge.nll2, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Brent",
            lower = c(logmean = 0.1),
            upper = c(logmean = 10))

fitmean = exp(coef(fit)[["logmean"]])
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa)*n,
                 fit_gamma = derlang2(time, D=fitmean)*n)
ggplot(df) + 
  geom_histogram(aes(x=interval)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Erlang after fit (2)"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> Erlang")


######################################################################
## gamma -> Pseudo Erlang
######################################################################

dperlang <- function(x, mean, kappa, shape=12, log=FALSE) { # shape=n; mean=D
  r <- kappa2r(kappa, shape)
  a <- (1-1/r^shape)/(mean*(1-1/r))
  density <- 0
  
  for (i in 1:shape) {
    innerp <- 1
    
    for (j in 1:shape) {
      if (j != i) {
        innerp <- innerp * r^(j-1)/(r^(j-1) - r^(i-1))
      }
    }
    density <- density + innerp * a * r^(i-1) * exp(-x*a*r^(i-1))
  }
  
  if (log) {
    density <- log(density)
  }
  return(density)
} 

gPE.nll <- function(logmean, logkappa, interval) {
  -sum(dperlang(interval, mean = exp(logmean), kappa = exp(logkappa), log=TRUE))
}

startPar = list(logmean = 2, logkappa = -1)
fit <- mle2(gPE.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = exp(coef(fit)[["logmean"]])
fitkappa = exp(coef(fit)[["logkappa"]])
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa)*n,
                 fit_gamma = dperlang(time, mean=fitmean, kappa=fitkappa)*n)
ggplot(df) + 
  geom_histogram(aes(x=interval)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Pseudo Erlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> PErlang")

## lager n -> (optim) function cannot be evaluated at initial parameters?
## shape difference between Gamma and Pseudo Erlang (same mean & kappa)
## PErlang fit ~= Erlang < Gamma
quit()
###############
kappa=0.3
mean=7
time <- timeSeq(1, 50, FALSE)
df <- data.frame(Time = time,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 perlang = dperlang(time, mean, kappa))
ggplot(df, aes(x=Time)) +
  geom_line(aes(y=gamma, color='Gamma'), linewidth=1) +
  geom_line(aes(y=perlang, color='Pseudo Erlang'), linewidth=1)

num_mean <- sum(df$Time * df$perlang * 1) # ts = 1
num_var <- sum(df$Time^2 * df$perlang * 1) - num_mean^2
num_kappa <- num_var / num_mean^2
num_kappa









