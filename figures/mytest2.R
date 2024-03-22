library(pseudoErlang)
library(deSolve)
library(bbmle)
library(ggplot2)

library(shellpipes)
startGraphics()

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

######################################################################
# hard test (large kappa) 

kappa <- 0.8 # maximum = 0.97
nfix <- 12

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
#print(fitW$fit) 
plotFit(fitW, df, fixedPar, type="SIgR", title = paste0("PE (n=12) -> PE (n=12) | kappa = " , kappa))


