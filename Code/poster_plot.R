# Poster Plots
library(deSolve)
library(ggplot2); theme_set(theme_classic(base_size = 23))
library(bbmle)
##library(devtools)
##install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(patchwork)

###################################################################### expDis

pdf("your_plot.pdf", width = 7, height = 5)
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

p <- ggplot(df, aes(x=Time, y = PD, color = n)) + 
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 1/γ, linetype = "dashed", linewidth = 1, color="black") +
  scale_color_manual(values = c('red', 'blue', 'green', 'yellow2', 'pink', 'magenta4')) + 
  labs(x="Stage duration, days", y = "Probability Density")
p
ggsave("/Users/ningruixie/Desktop/Uni/BIO_4C12/OE3C/LCT.pdf", plot = p, width = 9, height = 5, units = "in")

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

##
df <- data.frame(Time = time, LCT4 = dfn4$LCT, LCT8 = dfn8$LCT,
                 GCT4 = dfn4$GCT, GCT8 = dfn8$GCT)
p <- ggplot(df, aes(x=Time)) +
  geom_line(aes(y=LCT4, color = "n = 4"), linewidth=1.5) +
  geom_line(aes(y=LCT8, color = "n = 8"), linewidth=1.5) +
  scale_color_manual(values = c("n = 4" = "blue", "n = 8" = "orange")) +
  labs(color = "LCT", x = "Time, days",
       y = "Newly Infected Cases (%)") +
  geom_line(aes(y=GCT4, color = "black"), linewidth=1.5, linetype="dashed") +
  geom_line(aes(y=GCT8, color = "black"), linewidth=1.5, linetype="dashed")

ggsave("/Users/ningruixie/Desktop/Uni/BIO_4C12/OE3C/incComp.pdf", plot = p, width = 9, height = 5, units = "in")

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
p1 <- plotFit(fitW, df, fixedPar, title = "LCT (n=2) -> LCT (n=2)")
#print(fitW$fit)

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
p2 <- plotFit(fitW, df, fixedPar, type="SIgR", title = "GCT (n=2) -> GCT (n=2)")

P <- p1 | p2 
ggsave("/Users/ningruixie/Desktop/Uni/BIO_4C12/OE3C/sameTfit.pdf", plot = P, width = 15, height = 5, units = "in")

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
p <- ggplot(data2, aes(x=n, y = RMSLE*100)) + 
  geom_line(color='black') + #linetype='dashed'
  geom_point(shape=23, color="black", fill="#FFA500", size=3) + #69b3a2
  geom_hline(yintercept=data2$RMSLE[1]*100, color="#FFA500", linewidth=1) +
  geom_hline(yintercept=data$RMSLE[1]*100, color="#0077BB", linewidth=1) +
  geom_line(data=data, aes(n, y = RMSLE*100), color='black') +
  geom_point(data=data, aes(n, y = RMSLE*100), shape=21, color="black", fill="#0077BB", size=3) +
  labs(title = "LCT vs. GCT", y = "Root Mean Square Log Error (%)", x = "Substage Numbers")

ggsave("/Users/ningruixie/Desktop/Uni/BIO_4C12/OE3C/errorComp.pdf", plot = p, width = 8, height = 5, units = "in")


