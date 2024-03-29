library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 12))
library(bbmle)
library(pseudoErlang)
library(patchwork)

library(shellpipes)
startGraphics(height=3, width=9)

loadEnvironments()

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
  startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
  
  plots <- list()
  for (i in names(fitPar)) {
    data <- data.frame(variable = integer(),
                       estimate = numeric(),
                       lower = numeric(),
                       upper = numeric())
    for (n in cnfit) {
      fixedPar <- list(n = n, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
      fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
      
      log_est <- coef(fitW$fit)[[paste0("log", i)]]
      log_stde <- coef(summary(fitW$fit))[paste0("log", i), 'Std. Error']
      lower <- exp(log_est - 1.96*log_stde)
      upper <- exp(log_est + 1.96*log_stde)
      data <- rbind(data, data.frame(variable=n, estimate=exp(log_est), lower=lower, upper=upper))
    }
    
    p <- ggplot(data, aes(x=variable, y=estimate)) + 
      geom_point() +
      coord_flip()+
      geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
      geom_hline(yintercept=fitPar[[i]], col='red', linewidth=1) + 
      labs(x="Fitting Substages", y = paste0("Estimate ", ifelse(i == "β", "beta", i)))
    plots <- c(plots, list(p))
  }
  return(plots)
}

nfix <- 5
cnfit <- seq(5,10)
fitPar = c(β=0.2, D=10, kappa=2/9)
plots <- parEst(nfix, cnfit, fitPar)

######################################################################

plots[[1]] | plots[[2]] | plots[[3]]


