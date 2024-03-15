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
gE.nll <- function(logmean, logkappa, interval) { 
  ts=max(interval)/length(interval)
  sinr <- as.data.frame(SInRFlow(β=0, D=exp(logmean), n=round(1/exp(logkappa)), μ=0, 
                                 N=1, I0=1, 
                                 ts=ts, T=max(interval)))
  -sum(dnorm(x=interval, mean=diff(sinr$R)/ts*1000, log=TRUE))
}

#sinr = as.data.frame(SInRFlow(0,7,4,0,1,1,ts,25))
#plot(timeSeq(ts,25), diff(sinr$R)/ts*n)

startPar = list(logmean = 2, logkappa = -1)
fit <- mle2(gE.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

time = seq(0, max(g_interval), length.out=n)
fitmean = exp(coef(fit)[["logmean"]])
fitkappa = exp(coef(fit)[["logkappa"]])
ts = max(g_interval)/length(g_interval)
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa)*n,
                 fit_gamma = diff(SInRFlow(β=0, D=fitmean, n=round(1/fitkappa), μ=0, 
                                      N=1, I0=1, 
                                      ts=ts, T=max(g_interval))[,"R"])/ts*1000)
ggplot(df) + 
  geom_histogram(aes(x=interval)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Gamma after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> Erlang")

######################################################################
## gamma -> Pseudo Erlang
######################################################################
gPE.nll <- function(logmean, logkappa, interval) { 
  ts=max(interval)/length(interval)
  sigr <- as.data.frame(sinnerFlow(β=0, D=exp(logmean), kappa=exp(logkappa), n=7, μ=0, 
                                 N=n, I0=n, 
                                 ts=ts, T=max(interval)))
  -sum(dnorm(x=interval, mean=diff(sigr$R)/ts, log=TRUE))
}

sigr = sinnerFlow(0,10,0.3,7,0,1000,1000,0.1,25)
plot(x=timeSeq(0.1,25),y=diff(sigr[,"R"])/0.1,type='l')
#(β, D, kappa, n, μ, N, I0, ts, T)

startPar = list(logmean = 2, logkappa = -1)
fit <- mle2(gPE.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

time = seq(0, max(g_interval), length.out=n)
fitmean = exp(coef(fit)[["logmean"]])
fitkappa = exp(coef(fit)[["logkappa"]])
ts = max(g_interval)/length(g_interval)
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa)*n,
                 fit_gamma = diff(sinnerFlow(β=0, D=fitmean, kappa=fitkappa, n=7, μ=0, 
                                           N=n, I0=n, 
                                           ts=ts, T=max(g_interval))[,"R"])/ts)
ggplot(df) + 
  geom_histogram(aes(x=interval)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Gamma after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> Pseudo Erlang")








