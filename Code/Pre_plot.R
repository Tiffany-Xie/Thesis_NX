library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 18))
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
time <- timeSeq(ts, T)
sinner <- sinnerFlow(β, D, kappa, fixn, μ, N=N, I0=I0, ts, T)
sinr <- SInRFlow(β, D, n, μ, N=N, I0=I0, ts, T)
dfn4 <- data.frame(Time = time, 
                   LCT = diff(sinr[,"Cinc"])/ts,
                   GCT = diff(sinner[,"Cinc"])/ts)
n <- 8
kappa <- 1/n
time <- timeSeq(ts, T)
sinner <- sinnerFlow(β, D, kappa, fixn, μ, N=N, I0=I0, ts, T)
sinr <- SInRFlow(β, D, n, μ, N=N, I0=I0, ts, T)
dfn8 <- data.frame(Time = time, 
                   LCT = diff(sinr[,"Cinc"])/ts,
                   GCT = diff(sinner[,"Cinc"])/ts)

#####
df <- data.frame(Time = time, LCT4 = dfn4$LCT, LCT8 = dfn8$LCT,
                 GCT4 = dfn4$GCT, GCT8 = dfn8$GCT)
ggplot(df, aes(x=Time)) +
  geom_line(aes(y=LCT4, color = "n = 4"), linewidth=1) +
  geom_line(aes(y=LCT8, color = "n = 8"), linewidth=1) +
  scale_color_manual(values = c("n = 4" = "blue", "n = 8" = "orange")) +
  labs(color = "LCT", x = "Time, days",
       y = "Newly Infected Cases (proportion)") +
  geom_line(aes(y=GCT4, color = "black"), linewidth=1, linetype="dashed") +
  geom_line(aes(y=GCT8, color = "black"), linewidth=1, linetype="dashed")
#####
dfn4$id <- 'n=4'
dfn8$id <- 'n=8'
combined_df <- rbind(dfn4, dfn8)
long_df <- combined_df %>% 
  pivot_longer(cols = c(LCT, GCT), names_to = "variable", values_to = "value")


ggplot(long_df, aes(x = Time, y = value, color = id, linetype = variable)) +
  geom_line(linewidth=1) +
  scale_color_manual(values = c('n=4' = 'blue', 'n=8' = 'orange')) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(color = "n", linetype = "Model Type",
       x = "Time, days",
       y = "Newly Infected Cases (proportion)")

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
plotFit(fitW, df, fixedPar, title = "LCT (n=2) -> LCT (n=2)")
#print(fitW$fit)

# fit with different n
###################################################################### 

fixedPar <- list(n = 6, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) # what if: initial pop <<>> actual pop
fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotFit(fitW, df, fixedPar, title = "LCT (n=2) -> LCT (n=6)")

fixedPar <- list(n = 10, μ = μ, N = N, I0 = I0, ts = ts,  nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotFit(fitW, df, fixedPar, title = "E (n=2) -> E (n=10)")

fixedPar <- list(n = 12, μ = μ, N = N, I0 = I0, ts = ts,  nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll)
plotFit(fitW, df, fixedPar, title = "LCT (n=2) -> LCT (n=12)")

######################################################################
# Fit Pseudo Erlang with Pseudo Erlang
###################################################################### 

kappa <- 5/9
nfix <- 2

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)

nfit <- 2
startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
#print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "GCT (n=2) -> GCT (n=2)")

# should I be worried?
nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
#print(fitW$fit)
plotFit(fitW, df, fixedPar, type="SIgR", title = "GCT (n=2) -> GCT (n=12)")


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
plotFit(fitW, df, fixedPar, type="SIgR", title = "LCT (n=2) -> GCT (n=6)")

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "LCT (n=2) -> GCT (n=12)")

######################################################################
# show that the observed data, when changed, can also be fitted using nfix=12

n <- 7
sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)
df = simObs(sinr, arp, nbs, seed = 1001)

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
plotFit(fitW, df, fixedPar, type="SIgR", title = "LCT (n=7) -> GCT (n=12)")

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

## save time version...

β <- 0.2
D <- 10
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200
nbs<-1000
arp<-0.9
kappa=2/9

######################################################################

parEst <- function(nfix, cnfit, fitPar) {
  sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
  df <- simObs(sigr, arp, nbs, seed = 73)
  startPar_g <- list(logβ=-1, logD=2, logkappa=-0.1)
  startPar_l <- list(logβ=-1, logD=2)
  fixedPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T)
  
  for (par in names(fitPar)) {
    data <- data.frame(variable = integer(),
                       estimate_g = numeric(),
                       lower_g = numeric(),
                       upper_g = numeric(),
                       estimate_l = numeric(),
                       lower_l = numeric(),
                       upper_l = numeric())
    assign(paste("data", par, sep="_"), data)
  }
  
  for (n in cnfit) {
    fitW_g <- simplFit(startPar_g, fixedPar, df, sir.nll.g)
    fitW_l <- simplFit(startPar_l, fixedPar, df, sir.nll)
    
    log_est_beta_l <- coef(fitW_l$fit)[["logβ"]]
    log_est_beta_g <- coef(fitW_g$fit)[["logβ"]]
    log_stde_beta_l <- coef(summary(fitW_l$fit))[paste0("log", i), 'Std. Error']
    lower_beta_l <- exp(log_est - 1.96*log_stde)
    upper <- exp(log_est + 1.96*log_stde)
  }
}

parEst <- function(nfix, cnfit, fitPar) {
  sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
  df <- simObs(sigr, arp, nbs, seed = 73)
  startPar_g <- list(logβ=-1, logD=2, logkappa=-0.1)
  startPar_l <- list(logβ=-1, logD=2)
  
  plots <- list()
  for (i in names(fitPar)) {
    data_g <- data.frame(variable = integer(),
                       estimate_g = numeric(),
                       lower_g = numeric(),
                       upper_g = numeric())
    data_l <- data.frame(variable = integer(),
                         estimate_l = numeric(),
                         lower_l = numeric(),
                         upper_l = numeric())
    for (n in cnfit) {
      fixedPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
      fitW_g <- simplFit(startPar_g, fixedPar, df, sir.nll.g)
      
      log_est <- coef(fitW_g$fit)[[paste0("log", i)]]
      log_stde <- coef(summary(fitW_g$fit))[paste0("log", i), 'Std. Error']
      lower <- exp(log_est - 1.96*log_stde)
      upper <- exp(log_est + 1.96*log_stde)
      data_g <- rbind(data_g, data.frame(variable=n, estimate=exp(log_est), lower=lower, upper=upper))
      
      if (i != "kappa") {
        fixedPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
        fitW_l <- simplFit(startPar_l, fixedPar, df, sir.nll)
        
        log_est <- coef(fitW_l$fit)[[paste0("log", i)]]
        log_stde <- coef(summary(fitW_l$fit))[paste0("log", i), 'Std. Error']
        lower <- exp(log_est - 1.96*log_stde)
        upper <- exp(log_est + 1.96*log_stde)
        data_l <- rbind(data_l, data.frame(variable=n, estimate=exp(log_est), lower=lower, upper=upper)) 
      }
    }
    print(data_g)
    print(data_l)
    #p <- ggplot(data, aes(x=variable, y=estimate)) + 
    #  geom_point() +
    #  coord_flip()+
    #  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
    #  geom_hline(yintercept=fitPar[[i]], col='red', linewidth=1) + 
    #  labs(x="Fitting Substages", y = paste0("Estimate ", ifelse(i == "β", "beta", i)))
    #plots <- c(plots, list(p))
  }
  #return(plots)
}

nfix <- 5
cnfit <- seq(5,10)
fitPar = c(beta=0.2, D=10, kappa=2/9)
plots <- parEst(nfix, cnfit, fitPar)

######################################################################
# new MSE plot ??
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


multRMSLE <- function(data2fit, startPar,type="SInR", n2fit) {
  rmsledf <- data.frame(n = n2fit)
  RMSLE <- c()
  for (n in n2fit) {
    if (type=="SInR") {
      fPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T)
      fitW <- simplFit(startPar, fPar, data2fit, sir.nll)
      logβ <- coef(fitW$fit)[["logβ"]]
      logD <- coef(fitW$fit)[["logD"]]
      mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(logβ), D=exp(logD), 
                                                       n=fPar$n, μ=fPar$μ, N=fPar$N, 
                                                       I0=fPar$I0, ts=fPar$ts, T=fPar$T)))
    }
    else if (type=="SIgR") {
      fPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T)
      fitW <- simplFit(startPar, fPar, df, sir.nll.g)
      logβ <- coef(fitW$fit)[["logβ"]]
      logD <- coef(fitW$fit)[["logD"]]
      logkappa <- coef(fitW$fit)[["logkappa"]]
      mod.prep <- as.data.frame(as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa),
                                                         n=fPar$n, μ=fPar$μ, N=fPar$N, 
                                                         I0=fPar$I0, ts=fPar$ts, T=fPar$T)))
    }
    fitInc = diff(mod.prep$Cinc)
    #mse <- mean((data2fit$inc - fitInc)^2)
    rmsle <- sqrt(sum(log((fitInc+1)/(df$inc+1))^2)/dim(df)[1])
    RMSLE <- c(RMSLE, rmsle)
  }
  rmsledf["RMSLE"] = RMSLE
  return(rmsledf)
}

startPar <- list(logβ=-1, logD=2)
n2fit <- c(2,4,6,8,10)
data <- multRMSLE(df, startPar,type="SInR", n2fit) 

ggplot(data, aes(x=n, y = RMSLE*100)) + 
  geom_hline(yintercept=data$RMSLE[1]*100, color="#0077BB", linewidth=1) +
  geom_line(color='black') +
  geom_point(shape=21, color="black", fill="#0077BB", size=3) +
  labs(title = "LCT", y = "Root Mean Square Log Error (%)", x = "Substage Numbers")


#####
kappa <- 5/9
nfix <- 2

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)


startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
data2 <- multRMSLE(df, startPar,type="SIgR", n2fit) 
ggplot(data2, aes(x=n, y = RMSLE*100)) + 
  geom_line(color='black') + #linetype='dashed'
  geom_point(shape=23, color="black", fill="#FFA500", size=3) + #69b3a2
  geom_hline(yintercept=data2$RMSLE[1]*100, color="#FFA500", linewidth=1) +
  geom_hline(yintercept=data$RMSLE[1]*100, color="#0077BB", linewidth=1) +
  geom_line(data=data, aes(n, y = RMSLE*100), color='black') +
  geom_point(data=data, aes(n, y = RMSLE*100), shape=21, color="black", fill="#0077BB", size=3) +
  labs(title = "LCT vs. GCT", y = "Root Mean Square Log Error (%)", x = "Substage Numbers")







