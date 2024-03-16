library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 17))
library(bbmle)
##library(devtools)
##install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)

######################################################################

erlang <- function(x, n, γ) {
  (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
} 

time <- timeSeq(0.1, 40)
N <- c(1,2,4,10,20,70)
γ <- 0.1

df <- as.data.frame(expand.grid(Time=time,n=paste0("n=",N)))
PD <- c()
for (n in N) {
  PD <- c(PD, erlang(time, n, γ))
}
df$PD <- PD

ggplot(df, aes(x=Time, y = PD, color = n)) + 
  geom_line(linewidth = 1) +
  scale_color_manual(values = c('red', 'blue', 'green', 'yellow2', 'pink', 'magenta4')) + 
  labs(x="Stage duration, days", y = "Probability Density")
  



######################################################################


compareSDD <- function(D, n, nfix, ts, T, otherPar) {
  sinr <- SInRFlow(otherPar[["β"]], D, n, 
                   otherPar[["μ"]], otherPar[["N"]], 
                   otherPar[["I0"]], ts, T)
  sinner <- sinnerFlow(otherPar[["β"]], D, kappa=1/n, nfix, 
                       otherPar[["μ"]], otherPar[["N"]],
                       otherPar[["I0"]], ts, T)
  
  df <- data.frame(time = timeSeq(ts, T),
                   SInR = diff(sinr[,"R"])/ts,
                   SIgR = diff(sinner[,"R"])/ts)
  #flow =diff(sinr[,"R"])
  #mu <- sum(flow*df$time)
  #S <- sum(flow*(df$time)^2)
  #numkappa <- S/mu^2 - 1
  
  ggplot(df, aes(x = time)) + 
    geom_line(aes(y = SInR, color = "SInR"), linewidth = 1) + 
    geom_line(aes(y = SIgR, color = "SIgR"), linewidth = 1) + 
    labs(title = paste0("n =", n, " | nfix=", nfix), #, " | numD=", round(mu, 3), " | numK=", round(numkappa, 3)
         x = "Stage Duration, (days)",
         y = "Probability Density")
}

D = 10; nfix = 12; ts = 0.1; T = 40
otherPar = c(β = 0, μ = 0, N = 1, I0 = 1)
n = c(2,4,8,10)

plots <- lapply(n, function(n) compareSDD(D, n, nfix, ts, T, otherPar))
grid.arrange(grobs = plots, ncol = 2, nrow = 2)

######################################################################  Versatile_n

nn <- c(2,4,8)
D <- 10
ts <- 0.1
T <- 40
β = 0; μ = 0
N = 1; I0 = 1
time <- timeSeq(ts, T)

df <- expand.grid(Time = time, n = nn)

Sdd <- c()
for (n in nn) {
  sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
  sdd <- diff(sinr[,"R"])/ts
  Sdd <- c(Sdd, sdd)
}

df[["Sdd"]] <- Sdd

ggplot(df, aes(x=Time, y=Sdd, color = factor(n))) + geom_line(linewidth = 1) +
  geom_vline(xintercept = mu, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Stage duration (days)",
       y = "Probability Density")


######################################################################
library(tidyr)
library(dplyr)
β <- 0.2
D <- 10
fixn <- 12
μ <- 0.001
ts <- 0.1
T <- 100

N = 1
I0 = 0.01

n <- 4
kappa <- 1/n
time <- timeSeq(ts, T, FALSE)
sinner <- sinnerFlow(β, D, kappa, fixn, μ, N=N, I0=I0, ts, T)
sinr <- SInRFlow(β, D, n, μ, N=N, I0=I0, ts, T)
dfn4 <- data.frame(Time = time, 
                   Flat = rowSums(sinr[,3:(ncol(sinr)-2)]),
                   Geom = rowSums(sinner[,3:(ncol(sinner)-2)]))
n <- 8
kappa <- 1/n
time <- timeSeq(ts, T, FALSE)
sinner <- sinnerFlow(β, D, kappa, fixn, μ, N=N, I0=I0, ts, T)
sinr <- SInRFlow(β, D, n, μ, N=N, I0=I0, ts, T)
dfn8 <- data.frame(Time = time, 
                   Flat = rowSums(sinr[,3:(ncol(sinr)-2)]),
                   Geom = rowSums(sinner[,3:(ncol(sinner)-2)]))

dfn4$id <- 'n=4'
dfn8$id <- 'n=8'
combined_df <- rbind(dfn4, dfn8)
long_df <- combined_df %>% 
  pivot_longer(cols = c(Flat, Geom), names_to = "variable", values_to = "value")


ggplot(long_df, aes(x = Time, y = value, color = id, linetype = variable)) +
  geom_line(linewidth=1) +
  scale_color_manual(values = c('n=4' = 'blue', 'n=8' = 'orange')) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(color = "n", linetype = "Model Type",
       x = "Time, (days)",
       y = "Population Density",
       title = "Infctious Population Dynamics")

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
# run tempfunc before

sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 73)

startPar <- list(logβ=-1, logD=2)
fixedPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T)
# Consider estimating I0 (while fixing total pop size)

fitW <- simplFit(startPar, fixedPar, df, sir.nll)
#plotTrace(fitW)
plotFit(fitW, df, fixedPar, title = "SInR (n=2) -> SInR (n=2)")
#print(fitW$fit)

# fit with different n
###################################################################### 

fixedPar <- list(n = 6, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) # what if: initial pop <<>> actual pop
fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotFit(fitW, df, fixedPar, title = "SInR (n=2) -> SInR (n=6)")

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
#print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "SIgR (n=12) -> SIgR (n=12)")

# should I be worried?
nfit <- 6
startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
#print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "SIgR (n=12) -> SIgR (n=6)")


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
plotFit(fitW, df, fixedPar, type="SIgR", title = "SInR (n=2) -> SIgR (n=6)")

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "SInR (n=2) -> SIgR (n=12)")

######################################################################
# show that the observed data, when changed, can also be fitted using nfix=12

n <- 7
sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 1001)

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "SInR (n=7) -> SIgR (n=12)")

######################################################################
β <- 0.2; D <- 10; μ <- 0.01  
N <- 10000; I0 <- 10
ts <- 1; T <- 200
df <- data.frame(Time = timeSeq(ts, T))

n <- 2 
sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df[["n_2"]] = diff(sinr[,"Cinc"])

n <- 7
sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df[["n_7"]] = diff(sinr[,"Cinc"])

ggplot(df, aes(x=Time)) + 
  geom_line(aes(y=n_2), color = 'black', linewidth = 1) +
  geom_line(aes(y=n_7), color = "pink2", linewidth = 1) +
  labs(y = "Incidence", title = "Difference Between Incidence Data")
  #scale_color_manual(values = c('n_2' = 'blue', 'n_7' = 'orange'))

######################################################################
# error bar plot
######################################################################

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

###

data <- errplot(kappa=2/9, nfix=6, seq(2,10,by=2))

ggplot(data, aes(x=variable, y=estimate)) + 
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  geom_hline(yintercept=2/9, col='red', linewidth=1) + 
  labs(x="Fitting Substages", y = "Estimate κ") +
  annotate("text", x = 2, y = 2/9, label = "2/9", vjust = -1, color = "red", size=5)






