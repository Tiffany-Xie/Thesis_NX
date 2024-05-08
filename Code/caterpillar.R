#Caterpillar Plot
library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 18))
library(bbmle)
library(pseudoErlang)
library(metafor)
library(dplyr)
library(patchwork)

######################################################################

set.seed(5132)
k <- 250
vi <- c(0.2, 1) #rchisq(k, df=1) * .03
yi <- c(3,7) #rnorm(k, rnorm(k, 0.5, 0.4), sqrt(vi))

### fit RE model
res <- rma(yi, vi)

### create plot
forest(yi, vi,
       xlim=c(2,8),        ### adjust horizontal plot region limits
       order="obs",             ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       col="gray40",            ### change color of points/CIs
       psize=2,                 ### increase point size
       cex.lab=1, cex.axis=1,   ### increase size of x-axis title/labels
       lty=c("solid","blank"))  ### remove horizontal line at top of plot)
abline(v=4, lty=2)

### draw points one more time to make them easier to see
points(sort(yi), k:1, pch=19, cex=0.5)

### add summary polygon at bottom and text
addpoly(res, mlab="", cex=1)
text(-2, -2, "RE Model", pos=4, offset=0, cex=1)

######################################################################
### simulate data for the first set
set.seed(5132)
k <- 250
vi1 <- rchisq(k, df=1) * .03
yi1 <- rnorm(k, rnorm(k, 0.5, 0.4), sqrt(vi1))

### simulate data for the second set
set.seed(1234)  # set a different seed for reproducibility
vi2 <- rchisq(k, df=1) * .03
yi2 <- rnorm(k, rnorm(k, 0.7, 0.3), sqrt(vi2))

### combine the data
vi <- c(vi1, vi2)
yi <- c(yi1, yi2)

### fit RE model
res <- rma(yi, vi)

group <- rep(c("Group 1", "Group 2"), each = k)
### create plot
forest(yi, vi,
       xlim=c(-2.5,3.5),
       #order="obs",
       slab=group,  # use the group vector to label the studies
       annotate=TRUE,  # show study labels
       efac=0,
       pch=19,
       col=rep(c("red", "blue"), each = k),  # specify colors for each data point
       psize=2,
       cex.lab=1, cex.axis=1,
       lty=c("solid","blank"))

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



