library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 8))
library(bbmle)
library(pseudoErlang)
library(patchwork)

library(shellpipes)
startGraphics(height=5, width=9)

loadEnvironments()

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

sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 73)

######################################################################

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
org <- plotFit(fitW, df, fixedPar, type="SIgR", title = "LCT (n=2) -> GCT (n=12)")

######################################################################

nfit <- 6
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
cfitshape <- plotFit(fitW, df, fixedPar, type="SIgR", title = "LCT (n=2) -> GCT (n=6)")

######################################################################

n <- 7
sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 73)

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
cdata <- plotFit(fitW, df, fixedPar, type="SIgR", title = "LCT (n=7) -> GCT (n=12)")

######################################################################

(org | (cfitshape / cdata)) + plot_layout(widths = c(2.5, 1))






