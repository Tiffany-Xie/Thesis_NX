library(deSolve)
library(devtools)
#install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(ggplot2)

######################################################################
nPE <- 12
ts <- 0.1
T <- 200
mu <- 10
μ <- 0

###################################################################### GOOD

kappa <- 1/4 ## nE == 4
ode <- Integration(nPE, mu, kappa, ts, T)
gamm <- pgamma(timeSeq(ts, T, FALSE), 1/kappa, 1/(mu*kappa))

compPlot(gamm, ode, ts, T)
compPlotDens(gamm, ode, ts, T)

parCheck(gamm, ts, T)
parCheck(ode, ts, T)

###################################################################### BAD

parCheck_bad <- function(flow, ts, T) {
  flow <- diff(flow)
  time <- timeSeq(ts, T, FALSE)[-length(time)] # == +NA
  tot <- sum(flow)
  mu <- sum(flow*time)
  S <- sum(flow*time^2)
  kappa <- S/mu^2 - 1
  return(c(tot = tot, mu = mu, kappa = kappa))
}

parCheck_bad(gamm,ts,T)
parCheck_bad(ode,ts,T)
 
######################################################################

kappa <- 1/8 ## nE == 8
ode <- Integration(nPE, mu, kappa, ts, T)
gamm <- pgamma(timeSeq(ts, T, FALSE), 1/kappa, 1/(mu*kappa))

compPlot(gamm, ode, ts, T)
compPlotDens(gamm, ode, ts, T)

parCheck(gamm, ts, T)
parCheck(ode, ts, T)

######################################################################
kappa <- 1/2 ## nE == 2
ode <- Integration(nPE, mu, kappa, ts, T)
gamm <- pgamma(timeSeq(ts, T, FALSE), 1/kappa, 1/(mu*kappa))

compPlot(gamm, ode, ts, T)
compPlotDens(gamm, ode, ts, T)

parCheck(gamm, ts, T)
parCheck(ode, ts, T)



