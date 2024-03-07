library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)

##library(remotes)
##install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)

library(shellpipes)
## This statement will do nothing inside the pipeline, but allows to run this script interactively
rpcall("hardtest_error.Rout .pipestar Code/hardtest_error.R tempfunc.rda")

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


nfit <- 6
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
#print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = paste0("PE (n=12) -> PE (n=6) | kappa = " , kappa))

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
  startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
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

data <- errplot(kappa=2/9, nfix=12, c(5,6,7,8,9,10,11,12))

ggplot(data, aes(x=variable, y=estimate)) + 
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  geom_hline(yintercept=2/9, col='red', linewidth=1) + 
  labs(x="n fit", y = "Estimate Value", title="kappa=2/9, nfix=12")

###

data <- errplot(kappa=2/9, nfix=6, seq(2,10,by=2))

ggplot(data, aes(x=variable, y=estimate)) + 
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  geom_hline(yintercept=2/9, col='red', linewidth=1) + 
  labs(x="n fit", y = "Estimate Value", title="kappa=2/9, nfix=12")

###################################################################### 
# PE -> PE

kappa <- 2/9
nfix <- 12

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)

# should I be worried?
nfit <- 2
startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "PE (n=12) -> PE (n=2)")

nfit <- 3  
startPar <- list(logβ=-1, logD=2, logkappa=-1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "PE (n=12) -> PE (n=3)")


# Warning message:
# In mle2(likelihood.m, data = list(obs = obs), start = startPar,  :
#          couldn't invert Hessian



