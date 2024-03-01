library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
##library(devtools)
##install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)

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

startPar <- list(logβ=-1, logD=2)
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

# should I be worried?
nfit <- 2
startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "PE (n=12) -> PE (n=2)")

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
df = simObs(sinr, arp, nbs, seed = 73)

######################################################################
# show that changing the nfix won't affect fitting

nfit <- 6
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "E (n=2) -> PE (n=6)")

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "E (n=2) -> PE (n=12)")

######################################################################
# show that the observed data, when changed, can also be fitted using nfix=12

n <- 7
sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 1001)

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "E (n=7) -> PE (n=12)")


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





