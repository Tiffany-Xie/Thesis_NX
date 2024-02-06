library(deSolve)
library(devtools)
#install_github("Tiffany-Xie/pseudoErlang")
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

simplFit <- function(startPar, fixedPar, datadf, optMethod) {
  obs <- datadf$obs
  trace_betae <<- c()
  trace_De <<- c()
  
  fit0 <- mle2(sir.nll, 
               data=list(obs=obs),
               start=startPar, 
               fixed=fixedPar,
               method=optMethod,
               control=list(trace = 0))
  
   return(list(fit=fit0, trace_betae=trace_betae, trace_De=trace_De))
}

######################################################################

plotFit <- function(fitW, df, fPar) {
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
    labs(title = "Fitting Result", x = "Time, (days)", y = "Incidence")
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

######################################################################     

β <- 0.2
D <- 10
n <- 4
μ <- 0.01  
S0 <- 999
I0 <- 1
ts <- 1
T <- 100

arp <- 0.9
nbs <- 1000

###################################################################### 

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, seed=72)

startPar <- list(βe=-1, De=2)
fixedPar <- list(n = 1, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T) # what if: initial pop <<>> actual pop

fitW <- simplFit(startPar, fixedPar, df, "Nelder-Mead")
plotFit(fitW, df, fixedPar)
plotTrace(fitW)

## Change time ####################################################################

T <- 200

startPar <- list(βe=-0.5, De=2)
fixedPar <- list(n = n, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T)

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, seed=72)

fitW <- simplFit(startPar, fixedPar, df, "Nelder-Mead")
plotFit(fitW, df, fixedPar)
plotTrace(fitW)


## Change population (N) #################################################################### Change population (N)

T <- 100
S0 <- 99999

startPar <- list(βe=-1.2, De=2)
fixedPar <- list(n = n, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T)

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, seed=72)

fitW <- simplFit(startPar, fixedPar, df, "Nelder-Mead")
plotFit(fitW, df, fixedPar) # noise deviates from actual?
plotTrace(fitW)

## Change ts ####################################################################

S0 <- 999
ts <- 0.1 # works

startPar <- list(βe=-0.5, De=2)
fixedPar <- list(n = n, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T)

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, seed=72)

fitW <- simplFit(startPar, fixedPar, df, "Nelder-Mead")
plotFit(fitW, df, fixedPar) # noise deviates from actual?
plotTrace(fitW)


quit()
######################################################################









