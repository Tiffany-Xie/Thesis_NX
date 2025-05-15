library(deSolve)
library(ggplot2); theme_set(theme_classic(base_size = 13))
library(bbmle)
library(pseudoErlang)

library(shellpipes)
startGraphics(height=5, width=7)

loadEnvironments()

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


#####
kappa <- 5/9
nfix <- 2

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)


startPar <- list(logβ=-1, logD=2, logkappa=-0.5)
data2 <- multRMSLE(df, startPar,type="SIgR", n2fit) 
ggplot(data2, aes(x=n, y = RMSLE*100)) + 
  geom_line(color='black', linewidth=1.5) + #linetype='dashed'
  geom_hline(yintercept=data2$RMSLE[1]*100, color="red", linewidth=1.5) +
  geom_point(shape=23, color="black", fill="#FFA500", size=3) + #69b3a2
  #geom_hline(yintercept=data$RMSLE[1]*100, color="#0077BB", linewidth=1) +
  geom_line(data=data, aes(n, y = RMSLE*100), color='black', linewidth=1.5) +
  geom_point(data=data, aes(n, y = RMSLE*100), shape=21, color="black", fill="#0077BB", size=3) +
  labs(y = "Root Mean Square Log Error (%)", x = "Substage Numbers") +
  annotate("text", x = 9.7, y = 11.2, label = "LCT", color = "#0077BB", size = 7) +
  annotate("text", x = 9.7, y = 3, label = "GCT", color = "#FFA540", size = 7) +
  ylim(0, NA)
