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


nfit <- 2
startPar <- list(logβ=-1, logD=2, logkappa=-0.3)
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
# hard test (large kappa) 

kappa <- 0.97 # maximum = 0.97
nfix <- 12

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)

nfit <- 12  
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
print(fitW$fit) 
plotFit(fitW, df, fixedPar, type="SIgR", title = "PE (n=12) -> PE (n=12)")


nfit <- 6
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "PE (n=12) -> PE (n=6)")

######################################################################
# error bar plot

β <- 0.2
D <- 10
n <- 2
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200
nbs<-1000
arp<-0.9

errplot <- function(kappa, nfix, cnfit) { # cofInterv=0.95
  sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
  df <- simObs(sigr, arp, nbs, seed = 73)
  startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
  data <- data.frame(variable = integer(),
                   estimate = numeric(),
                   lower = numeric(),
                   upper = numeric())
  print(data)

  for (n in cnfit) {
      fixedPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
      fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
      
      log_est_k <- coef(fitW$fit)[["logkappa"]]
      log_stde_k <- coef(summary(fitW$fit))['logkappa', 'Std. Error']
      lower <- exp(log_est_k - 1.96*log_stde_k)
      upper <- exp(log_est_k + 1.96*log_stde_k)
      data <- rbind(data, data.frame(variable=n, estimate=exp(log_est_k), lower=lower, upper=upper))
      print(data)
    }
  return(data)
}  

data <- errplot(2/9, 12, c(5,6,7,8,9,10,11,12))

ggplot(data, aes(x=variable, y=estimate)) + 
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  geom_hline(yintercept=2/9, col='red', linewidth=1)



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





