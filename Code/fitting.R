library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)

## install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)

######################################################################

## Make a data frame with interval incidence, and then observations sampled from that
simObs <- function(sinr, arp, nbs, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  inc <- diff(sinr[,"Cinc"]) #/ts
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  df <- data.frame(Time = (sinr[,"time"][-1] + sinr[,"time"][-dim(sinr)[1]])/2,
                   inc = inc,
                   obs = obs)
  return(df)
}

######################################################################

## likelihood function with parameter-trajectory tracing 🙂
## FIXME: size needs to be _passed_ to this function √
sir.nll <- function(logβ, logD, n, μ, N, I0, ts, T, nbs, obs) {
  trace_logbeta <<- c(trace_logbeta, logβ)
  trace_logD <<- c(trace_logD, logD)
  
  out <- as.data.frame(SInRFlow(β=exp(logβ), D=exp(logD), n=n, μ=μ, N=N, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

## Similar function for geometric pseudo-Erlang (different params)
sir.nll.g <- function(logβ, logD, logkappa, n, μ, N, I0, ts, T, nbs, obs) {
  trace_logbeta <<- c(trace_logbeta, logβ)
  trace_logD <<- c(trace_logD, logD)
  trace_logkappa <<- c(trace_logkappa, logkappa)
  
  out <- as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa), n=n, μ=μ, N=N, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

## Wrapper function to fit the above likelihoods
simplFit <- function(startPar, fixedPar, datadf, likelihood.m, optMethod="Nelder-Mead") {
  obs <- round((datadf$obs)/arp)
  trace_logbeta <<- c()
  trace_logD <<- c()
  trace_logkappa <<- c()
  
  fit0 <- mle2(likelihood.m, 
               data=list(obs=obs),
               start=startPar, 
               fixed=fixedPar,
               method=optMethod,
               control=list(trace = 0))
  
  return(list(fit=fit0, trace_logbeta=trace_logbeta, trace_logD=trace_logD))
}

######################################################################

plotFit <- function(fitW, df, fPar, type="SInR", title = "Fitting Result") {
  logβ <- coef(fitW$fit)[["logβ"]]
  logD <- coef(fitW$fit)[["logD"]]
  
  if (type=="SInR") {
    mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(logβ), D=exp(logD), 
                                                     n=fPar$n, μ=fPar$μ, N=fPar$N, 
                                                     I0=fPar$I0, ts=fPar$ts, T=fPar$T)))
  }
  else if (type=="SIgR") {
    logkappa <- coef(fitW$fit)[["logkappa"]]
    mod.prep <- as.data.frame(as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa),
                                                     n=fPar$n, μ=fPar$μ, N=fPar$N, 
                                                     I0=fPar$I0, ts=fPar$ts, T=fPar$T)))
  }
  
  df["fitInc"] = diff(mod.prep$Cinc)
  mse <- mean((df$inc - df$fitInc)^2)
  
  p <- ggplot(df, aes(x=Time)) +
    geom_line(aes(y=obs, color='Observed')) +
    geom_line(aes(y=inc, color = 'Incidence'), linewidth=1) +
    geom_line(aes(y=fitInc, color = 'Fit Incidence'), linewidth=1.5, alpha=0.8) +
    labs(title = paste(title, "| MSE =", round(mse, digits=3)), x = "Time, (days)", y = "Incidence")
  return(p)
}


plotTrace <- function(fitW) {
  trace_logbeta <- fitW$trace_logbeta
  trace_logD <- fitW$trace_logD
  
  pardf <- data.frame(order = seq(0:(length(trace_logbeta)-1)), 
                      logbeta = trace_logbeta,
                      logD = trace_logD)
  
  ggplot(pardf, aes(x = logD, y = logbeta)) +
    geom_point(aes(color = order)) +
    geom_path(alpha = 0.5) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(x = "log(D)", y = "log(beta)", color = "Time")
}

######################################################################
# Fit Erlang with Erlang
######################################################################     

β <- 0.2
D <- 10
n <- 2
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200

arp <- 0.9
nbs <- 1000

###################################################################### 

sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 73)

## FIXME: don't call log β βe!! etc!! √
startPar <- list(logβ=-1.5, logD=2)
fixedPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T)
# Consider estimating I0 (while fixing total pop size)

fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotTrace(fitW)
plotFit(fitW, df, fixedPar, title = "E (n=2) -> E (n=2)")

print(fitW$fit)

# fit with different n
###################################################################### 

fixedPar <- list(n = 6, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) # what if: initial pop <<>> actual pop
fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotFit(fitW, df, fixedPar, title = "E (n=2) -> E (n=6)")

fixedPar <- list(n = 10, μ = μ, N = N, I0 = I0, ts = ts,  nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotFit(fitW, df, fixedPar, title = "E (n=2) -> E (n=10)")

fixedPar <- list(n = 12, μ = μ, N = N, I0 = I0, ts = ts,  nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotFit(fitW, df, fixedPar, title = "E (n=2) -> E (n=12)")

######################################################################
# Fit Pseudo Erlang with Pseudo Erlang
###################################################################### 

kappa <- 2/9
nfix <- 12

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)

nfit <- 12  
startPar <- list(logβ=-1, logD=2, logkappa=-1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "PE (n=12) -> PE (n=12)")


nfit <- 6 
startPar <- list(logβ=-1, logD=2, logkappa=-1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "PE (n=12) -> PE (n=6)")

######################################################################
# check shape difference (while maintaining mean) influence simulation

β <- 0.2
D <- 10
μ <- 0.01 
N <- 10000
I0 <- 10
ts <- 1
T <- 400

sinr2 <- SInRFlow(β, D, n=2, μ, N, I0, ts, T)
sinr8 <- SInRFlow(β, D, n=8, μ, N, I0, ts, T)
time <- timeSeq(ts, T)
plot(time, diff(sinr2[,"Cinc"]), 
     main = "Erlang (n=2,8): how changes in substages affect incidence", 
     type="l")
lines(time, diff(sinr8[,"Cinc"]), col="red")

######################################################################
# Fit Erlang with Pseudo Erlang
######################################################################
β <- 0.2
D <- 10
n <- 2
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200

sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 71)

######################################################################
# show that changing the nfix won't affect fitting

nfit <- 6
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "E (n=2) -> PE (n=6)")

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "E (n=2) -> PE (n=12)")

######################################################################
# show that the observed data, when changed, can also be fitted using nfix=12

n <- 10
sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 60)

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "E (n=10) -> PE (n=12)")



quit()
######################################################################

# change n doesn't seems to affect fitting a lot
# same for SIgR model
# change initial value will alter the fitting
# the effect of sub-stage number seems can be compensate by other parameter
# only observe somewhat large different when all the other parameters are same

# Erlang -> Pseudo Erlang did not meet expectations

# In mle2(likelihood.m, data = list(obs = obs), start = startPar,  :
# couldn't invert Hessian

# In mle2(likelihood.m, data = list(obs = obs), start = startPar,  :
# convergence failure: code=10 (degenerate Nelder-Mead simplex)






