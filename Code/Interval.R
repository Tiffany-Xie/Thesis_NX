library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)

set.seed(223)

numObs <- 1000
kappa <- 0.4
D <- 5

## k = 1/shape
## D = shape*scale = scale/kappa

v <- round(rgamma(numObs, shape=1/kappa, scale=D*kappa))

print(mean(v))
print(sd(v)^2/mean(v)^2)

df <- data.frame(v = v)
ggplot(df, aes(x = v)) + geom_histogram()

## start/end times

## start is a uniform
## end is start + v [above]
## round them both before subtracting

gammaDelay <- function(n, mean, kappa){
  
}

gammaDis <- function(time, mean, kappa) {
  alpha = 1/kappa
  beta = 1/(kappa*mean)
  beta^alpha / factorial(alpha-1) * time^(alpha-1) * exp(-beta*time)  
}


sir.nll <- function(logmean, logkappa, ts,T, obs) {
  #obs = round(obs)
  out <- gammaDis(timeSeq(ts, T, FALSE), exp(logmean), exp(logkappa))
  nll <- -sum(dnorm(x=obs, mean=out, sd=1, log=TRUE))
}

simplFit <- function(startPar, fixedPar, obs, likelihood.m, optMethod="Nelder-Mead") {
  
  fit0 <- mle2(likelihood.m, 
               data=list(obs=obs),
               start=startPar,
               fixed=fixedPar,
               method=optMethod,
               control=list(trace = 0))
  
  return(fit0)
}

ts = 1; T = 30; mean = 10; kappa = 1/2
g <- gammaDis(timeSeq(ts, T, FALSE), mean, kappa)

startPar = list(logmean = 2, logkappa = -0.5)
fixedPar = list(ts=1, T=30)
fit = simplFit(startPar, fixedPar, g, sir.nll)

df <- data.frame(time = timeSeq(ts, T, FALSE), 
                 obs=g, fitted = gammaDis(timeSeq(ts, T, FALSE), 
                                          exp(coef(fit)[["logmean"]]),
                                          exp(coef(fit)[["logkappa"]])))

ggplot(df, aes(time)) +
  geom_line(aes(y = obs, col = "obs")) + 
  geom_line(aes(y = fitted, col = "fitted")) + 
  labs(y = "Probability Density")


quit()
######################################################################

time <- 0:30
g <- gammaDis(time, 10, 0.3)
df <- data.frame(time=time, gamma = g)
ggplot(df, aes(x = time, y = gamma)) + geom_line()
  