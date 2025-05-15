#Caterpillar Plot
library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 18))
library(bbmle)
library(pseudoErlang)
library(metafor)
library(dplyr)
library(patchwork)

######################################################################

cVarEst <- function(sigr, sParL, sParG, fPar, seeds, par) {
  data_l <- list()
  for (seed in seeds) {
    df <- simObs(sigr, arp, nbs, seed = seed)
    fitWL <- simplFit(sParL, fPar, df, sir.nll)
    fitWG <- simplFit(sParG, fPar, df, sir.nll.g)
    parL <- coef(summary(fitWL$fit))
    krow <- matrix(data=c(log(1/fPar[["n"]]), 0, NA, NA), nrow=1,
                   dimnames = list("logkappa", c("Estimate", "Std. Error", "z value", "Pr(z)")))
    parL <- rbind(parL, krow)
    parG <- coef(summary(fitWG$fit))
    
    for (p in par) {
      parname <- paste0("log",p)
      estL <- parL[parname, "Estimate"]
      varL <- 1.96 * parL[parname, "Std. Error"]
      estG <- parG[parname, "Estimate"]
      varG <- 1.96 * parG[parname, "Std. Error"]
      data_l[[p]] <- bind_rows(data_l[[p]], data.frame(estL=exp(estL), varL=exp(estL)*(exp(varL)-1), estG=exp(estG), varG=exp(estG)*(exp(varG)-1)))
    }
    print(paste(seed,"Done"))
  }
  return(data_l)
}

######################################################################     

β <- 0.2
D <- 10
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200
nbs<-1000
arp<-0.9
kappa<- 0.3

nfix = 5
nfit = 8

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)

par <- c("β", "D", "kappa")
sParG <- list(logβ=-1, logD=2, logkappa=-1)
sParL <- list(logβ=-1, logD=2)
fPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T)
seeds <- seq(70,70+100)

data = cVarEst(sigr, sParL, sParG, fPar, seeds, par)

######################################################################     

catPlot <- function(data, par) {
  par(mfrow = c(1, 2))
  for (p in names(par)) {
    pdata <- data[[p]]
    forest(pdata$estL, pdata$varL,
           #xlim=c(par[[p]]-2,par[[p]]+2), 
           order="obs",             
           slab=NA, annotate=FALSE, 
           efac=0,                  
           pch=19,                  
           col="gray40",            
           psize=2,                 
           cex.lab=1, cex.axis=1,   
           lty=c("solid","blank"))  
    abline(v=par[[p]], lty=2)
    title(main = paste("LCT,", p))
    
    forest(pdata$estG, pdata$varG,
           #xlim=c(par[[p]]-2,par[[p]]+2), 
           order="obs",             
           slab=NA, annotate=FALSE, 
           efac=0,                  
           pch=19,                  
           col="gray40",            
           psize=2,                 
           cex.lab=1, cex.axis=1,   
           lty=c("solid","blank"))  
    abline(v=par[[p]], lty=2)
    title(main = paste("GCT,", p))
    #par(mfrow = c(1, 1))
  }
}

par <- c(β=0.2, D=10, kappa=0.3)

catPlot(data, par)

######################################################################     

linplot <- function(data, par) {
  for (p in names(par)) {
    print(ggplot(data[[p]], aes(x=1:dim(temp)[1])) +
            geom_line(aes(y=sort(estL), color="LCT"), linewidth=1) +
            geom_line(aes(y=sort(estG), color="GCT"), linewidth=1) +
            geom_hline(yintercept = par[[p]], color='red', linewidth=1.5) +
            labs(title=p, x="", y=paste("Estimate", p)))
  }
}

linplot(data, par)



