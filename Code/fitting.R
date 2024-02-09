library(deSolve)
library(devtools)
install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)

######################################################################

simObs <- function(sinr, arp, nbs, seed) {
  set.seed(seed)
  inc <- diff(sinr[,"inc"]) #/ts
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  df <- data.frame(Time = (sinr[,"time"][-1] + sinr[,"time"][-dim(sinr)[1]])/2,
                   inc = inc,
                   obs = obs)
  return(df)
}

######################################################################

sir.nll <- function(βe, De, n, μ, S0, I0, ts, T, obs) {
  
  trace_betae <<- c(trace_betae, βe)
  trace_De <<- c(trace_De, De)
  
  out <- as.data.frame(SInRFlow(β=exp(βe), D=exp(De), n=n, μ=μ, S0=S0, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=1000, log=TRUE))
}

sir.nll.g <- function(βe, De, kappae, n, μ, S0, I0, ts, T, obs) {
  
  trace_betae <<- c(trace_betae, βe)
  trace_De <<- c(trace_De, De)
  trace_kappae <<- c(trace_kappae, kappae)
  
  out <- as.data.frame(sinnerFlow(β=exp(βe), D=exp(De), kappa=exp(kappae), n=n, μ=μ, S0=S0, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=1000, log=TRUE))
}

simplFit <- function(startPar, fixedPar, datadf, likelihood.m, optMethod) {
  obs <- round((datadf$obs)/arp)
  trace_betae <<- c()
  trace_De <<- c()
  trace_kappae <<- c()
  
  fit0 <- mle2(likelihood.m, 
               data=list(obs=obs),
               start=startPar, 
               fixed=fixedPar,
               method=optMethod,
               control=list(trace = 0))
  
  return(list(fit=fit0, trace_betae=trace_betae, trace_De=trace_De))
}

######################################################################

plotFit <- function(fitW, df, fPar, title = "Fitting Result") {
  βe <- coef(fitW$fit)[["βe"]]
  De <- coef(fitW$fit)[["De"]]
  
  mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(βe), D=exp(De), 
                                                   n=fPar$n, μ=fPar$μ, S0=fPar$S0, 
                                                   I0=fPar$I0, ts=fPar$ts, T=fPar$T)))
  df["fitInc"] = diff(mod.prep$inc)
  
  p <- ggplot(df, aes(x=Time)) +
    geom_line(aes(y=obs, color='Observed')) +
    geom_line(aes(y=inc, color = 'Incidence'), linewidth=1) +
    geom_line(aes(y=fitInc, color = 'Fit Incidence'), linewidth=1.5, alpha=0.8) +
    labs(title = title, x = "Time, (days)", y = "Incidence")
  return(p)
}


plotTrace <- function(fitW) {
  trace_betae <- fitW$trace_betae
  trace_De <- fitW$trace_De
  
  pardf <- data.frame(order = seq(0:(length(trace_betae)-1)), 
                      betae = trace_betae,
                      De = trace_De)
  
  ggplot(pardf, aes(x = De, y = betae)) +
    geom_point(aes(color = order)) +
    geom_path(alpha = 0.5) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(x = "log(D)", y = "log(beta)", color = "Time")
}

# Fit Erlang with Erlang
######################################################################
######################################################################     

β <- 0.2
D <- 10
n <- 4
μ <- 0.01  
S0 <- 999
I0 <- 1
ts <- 1
T <- 200

arp <- 0.9
nbs <- 1000

###################################################################### 

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 73)

startPar <- list(βe=-1.5, De=2)
fixedPar <- list(n = 4, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T) # what if: initial pop <<>> actual pop

fitW <- simplFit(startPar, fixedPar, df, sir.nll, "Nelder-Mead")
plotTrace(fitW)
plotFit(fitW, df, fixedPar, title = "Fitting Result (n = 4)")

# with different n
###################################################################### 

fixedPar <- list(n = 6, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T) # what if: initial pop <<>> actual pop
fitW <- simplFit(startPar, fixedPar, df, sir.nll, "Nelder-Mead")
plotFit(fitW, df, fixedPar, title = "SInR Fitting Result (n = 6)")

fixedPar <- list(n = 10, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T) # what if: initial pop <<>> actual pop
fitW <- simplFit(startPar, fixedPar, df, sir.nll, "Nelder-Mead")
plotFit(fitW, df, fixedPar, title = "SInR Fitting Result (n = 10)")

fixedPar <- list(n = 12, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T) # what if: initial pop <<>> actual pop
fitW <- simplFit(startPar, fixedPar, df, sir.nll, "Nelder-Mead")
plotFit(fitW, df, fixedPar, title = "SInR Fitting Result (n = 12)")


# Fit Pseudo Erlang with Pseudo Erlang
######################################################################
###################################################################### 

tempFit <- function(sigr, df, startP, fixedP, plot=TRUE) {
  fitW <- simplFit(startP, fixedP, df, sir.nll.g, "Nelder-Mead")
  
  βe <- coef(fitW$fit)["βe"]
  De <- coef(fitW$fit)["De"]
  kappae <- coef(fitW$fit)["kappae"]
  mod.prep <- as.data.frame(as.data.frame(sinnerFlow(β=exp(βe), D=exp(De), kappa=exp(kappae),
                                                     n=fixedP$n, μ=fixedP$μ, S0=fixedP$S0, 
                                                     I0=fixedP$I0, ts=fixedP$ts, T=fixedP$T)))
  df["fitInc"] = diff(mod.prep$inc)
  
  if (plot) {
    print(ggplot(df, aes(x=Time)) +
    geom_line(aes(y=obs, color='Observed')) +
    geom_line(aes(y=inc, color = 'Incidence'), linewidth=1) +
    geom_line(aes(y=fitInc, color = 'Fit Incidence'), linewidth=1.5, alpha=0.5) +
    labs(title = paste0("SIgR, nfit =", fixedPar$n), x = "Time, (days)", y = "Incidence")
  )}
  
  return(fitW)
}

###################################################################### 

kappa <- 2/9
nfix <- 12

sigr <- sinnerFlow(β, D, kappa, nfix, μ, S0, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)

nfit <- 12
startPar <- list(βe=-1, De=2, kappae=-1)
fixedPar <- list(n = nfit, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T) 

fitW <- tempFit(sigr, df, startPar, fixedPar,)

###################################################################### 
# check what's wrong
sigr12 <- sinnerFlow(β=exp(-1.619631), D=exp(2.323400), kappa=exp(-1.507136), n=12, μ, S0, I0, ts, T)
sigr6 <- sinnerFlow(β=exp(-1.614899), D=exp(2.319774), kappa=exp(-1.478050), n=6, μ, S0, I0, ts, T)

time <- timeSeq(ts, T)

plot(time, diff(sigr12[,"inc"]), main = "Pseudo Erlang (n=12,6), with slightly different parameters (post fit)")
lines(time, diff(sigr6[,"inc"]), type='l', col="red", lwd = 3)

######################################################################
# check shape difference (while maintaining mean) influence simulation
β <- 0.2
D <- 10
n <- 4
μ <- 0.01  
S0 <- 999
I0 <- 1
ts <- 1
T <- 400

sinr4 <- SInRFlow(β, D, n=4, μ, S0, I0, ts, T)
sinr8 <- SInRFlow(β, D, n=8, μ, S0, I0, ts, T)
time <- timeSeq(ts, T)
plot(time, diff(sinr4[,"inc"]), 
     main = "Erlang (n=4,8): how changes in substages affect incidence", 
     type="l")
lines(time, diff(sinr8[,"inc"]), col="red")

# Fit Erlang with Pseudo Erlang
######################################################################
######################################################################
β <- 0.2
D <- 10
n <- 4
μ <- 0.01  
S0 <- 999
I0 <- 1
ts <- 1
T <- 200

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 71)

nfit <- 12
startPar <- list(βe=-1, De=2, kappae=-0.1)
fixedPar <- list(n = nfit, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T) 

fitW <- simplFit(startPar, fixedPar, df, sir.nll.g, "Nelder-Mead")
plotTrace(fitW)
plotFit(fitW, df, fixedPar, title = "Fitting Erlang (n=4) with Pseudo Erlang (n=12)")



quit()
######################################################################

# change n doesn't seems to affect fitting a lot
# same for SIgR model
# change initial value will alter the fitting
# the effect of sub-stage number seems can be compensate by other parameter
# only observe somewhat large different when all the other parameters are same

# Erlang -> Pseudo Erlang did not meet expectations

# check & change titles







