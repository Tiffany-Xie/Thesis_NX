library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)

library(shellpipes)
loadEnvironments()
# source("./tempfunc.R")
######################################################################

## Continuous data should be fitted with derlang and dperlang
gammaDelay <- function(n, mean, kappa){
  v <- rgamma(n, shape=1/kappa, scale=mean*kappa)
  print(c(mean = mean, kappa = kappa,
          Dmean = mean(v), Dkappa = sd(v)^2/mean(v)^2))
  return(v)
}

## Discrete data should be fitted another way (make boxed likelihood functions)
gammaDelayDiscrete <- function(n, mean, kappa){
  v <- round(rgamma(n, shape=1/kappa, scale=mean*kappa))
  print(c(mean = mean, kappa = kappa,
          Dmean = mean(v), Dkappa = sd(v)^2/mean(v)^2))
  return(v)
}

## Boxed likelihood example (not working yet!!!)
blgamma <- function(n){
	## Do we need a special case for n=0?
	## should we check that it's an integer and/or round?
	return(pgamma(n+1/2, log=TRUE) - pgamma(n-1/2, log=TRUE))
}

######################################################################
## Gamma -> Gamma
######################################################################

## k = 1/shape
## D = shape*scale = scale/kappa

## start/end times
## start is a uniform
## end is start + v [above]
## round them both before subtracting

set.seed(223)
mean = 7
kappa = 0.1
n = 2000

g_interval = gammaDelay(n,mean,kappa) # v
startPar = list(logmean = 2.5, logkappa = -0.5)
time = seq(0, max(g_interval), length.out=n)

fit <- mle2(gg.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 fit_gamma = dgamma(time, shape=1/fitkappa, scale=fitmean*fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = ..density..)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Gamma after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> gamma")

######################################################################
## Gamma -> Erlang
######################################################################

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
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 fit_gamma = derlang(time, shape=round(1/fitkappa), rate=1/fitmean*round(1/fitkappa)))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = ..density..)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Erlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> Erlang (shape parameter included)")

######################################################################
## Gamma -> Pseudo Erlang
######################################################################

startPar = list(logmean = 2, logkappa = -1)
fit <- mle2(gPE.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 fit_gamma = dperlang(time, mean=fitmean, kappa=fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = ..density..)) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Pseudo Erlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "gamma -> PErlang")

######################################################################
## Random Pseudo Erlang numbers
######################################################################

set.seed(1001)
n <- 1000
mean <- 7
kappa <- 0.3

pe_interval <- rperlang(n, mean, kappa) # <<slow if c->inf>>
time = seq(0, max(pe_interval), length.out=n)
df <- data.frame(Time=time, interval=pe_interval,
                 perlang=dperlang(time, mean, kappa),
                 gamma=dgamma(time, 1/kappa, scale=mean*kappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = ..density..)) +
  geom_line(aes(x=Time, y=perlang, color="PErlang"), linewidth=1.5) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  labs(title="Random PErlang")

## sometimes higher than curve?
## count/curve inconsistency when using >><< kappa (similar with rgamma)
## best way to find the proposal distribution

######################################################################
## Pseudo Erlang -> Pseudo Erlang
######################################################################

startPar = list(logmean = 2.5, logkappa = -0.5)

fit <- mle2(gPE.nll, 
            data = list(interval = pe_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = pe_interval,
                 perlang = dperlang(time, mean, kappa),
                 fit_perlang = dperlang(time, fitmean, fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = ..density..)) +
  geom_line(aes(x=Time, y=perlang, color="PErlang"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_perlang, color="PErlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "PErlang -> PErlang")

######################################################################
## Pseudo Erlang -> Gamma
######################################################################

startPar = list(logmean = 2.5, logkappa = -0.5)

fit <- mle2(gg.nll, 
            data = list(interval = pe_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = pe_interval,
                 perlang = dperlang(time, mean, kappa),
                 fit_gamma = dgamma(time, 1/fitkappa, scale=fitmean*fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = ..density..)) +
  geom_line(aes(x=Time, y=perlang, color="PErlang"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Gamma after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = "PErlang -> Gamma")

######################################################################
## Test Pseudo Erlang CDF
######################################################################

mean <- 7
kappa <- 0.3
ts <- 0.1

time <- timeSeq(ts, 30)
perlangPDF <- dperlang(time, mean, kappa)
perlangCDF <- pperlang(time, mean, kappa)
numCDF <- numeric(length(time))
for (i in 1:length(time)) {
  numCDF[i] <- sum(ts * perlangPDF[1:i])
}

df <- data.frame(Time=time, numCDF=numCDF, aCDF=perlangCDF)
ggplot(df, aes(x=Time)) +
  geom_line(aes(y=numCDF, color='numeric CDF'), linewidth=1) +
  geom_line(aes(y=aCDF, color='actual CDF'), linewidth=1.5, linetype='dashed') +
  labs(title='Pseudo Erlang actual vs. numeric CDF', y="CD")


############### check
kappa=0.3
mean=7
time <- timeSeq(1, 1000, FALSE)
df <- data.frame(Time = time,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 perlang = dperlang(time, mean, kappa))
ggplot(df, aes(x=Time)) +
  geom_line(aes(y=gamma, color='Gamma'), linewidth=1.5, linetype='dashed') +
  geom_line(aes(y=perlang, color='Pseudo Erlang'), linewidth=1) +
  labs(title="Check: same mean & kappa")

ggplot(df[-1,], aes(x=Time)) +
  geom_line(aes(y=log(gamma), color='Gamma(l)'), linewidth=1.5, linetype='dashed') +
  geom_line(aes(y=log(perlang), color='Pseudo Erlang(l)'), linewidth=1) +
  labs(title="Check: same mean & kappa")

## area under the curve
sum(df$gamma*0.005)
sum(df$perlang*0.005)

quit()
num_mean <- sum(df$Time * df$perlang * 1) # ts = 1
num_var <- sum(df$Time^2 * df$perlang * 1) - num_mean^2
num_kappa <- num_var / num_mean^2
num_kappa

temp <- dperlang(time, 7, 0.02)
temp2 <- dgamma(time, shape=1/0.02, scale=7*0.02)

temp_df <- data.frame(Time=time, gamma=temp2, perlang=temp)
ggplot(temp_df, aes(x=Time)) +
  geom_line(aes(y=gamma,color='gamma')) +
  geom_line(aes(y=perlang,color='perlang'))

tempmean <- sum(temp2 * time * 0.012776)
tempvar <- sum(time^2 * temp2 * 0.012776) - tempmean^2
tempkappa <- tempvar/tempmean^2
tempkappa








