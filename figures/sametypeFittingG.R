library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 12))
library(bbmle)
library(pseudoErlang)
library(patchwork)

library(shellpipes)
startGraphics(height=7, width=9)

loadEnvironments()

######################################################################

onetimeParE <- function(fit, fitPar) { # fitPar = c(N = 1, D = 10)
  dfs <- list()
  plots <- list()
  
  for (i in names(fitPar)) {
    
    log_est <- coef(fitW$fit)[[paste0("log", i)]]
    log_stde <- coef(summary(fitW$fit))[paste0("log", i), 'Std. Error']
    lower <- exp(log_est - 1.96*log_stde)
    upper <- exp(log_est + 1.96*log_stde)
    
    df <- data.frame(variable = ifelse(i == "β", "beta", i),
                     estimate = exp(log_est),
                     lower = lower,
                     upper = upper)
    #print(df)
    
    #dfs[[paste("df", i, sep = "_")]] <- df
    
    p <- ggplot(df, aes(x=variable, y=estimate)) +
      geom_point() +
      geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
      geom_hline(yintercept=fitPar[[i]], col='red', linewidth=1) +
      labs(x = element_blank(), y = "Estimate Parameter")
    #print(p)
    
    plots <- c(plots, list(p))
  }
  #grid.arrange(grobs = plots, ncol = length(fitPar))
  return(plots)
}

######################################################################

β <- 0.2
D <- 10
nfix <- 12
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200
kappa <- 2/9

arp <- 0.9
nbs <- 1000

startPar <- list(logβ=-1, logD=2, logkappa=-1)

###################################################################### Same n

sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
df <- simObs(sigr, arp, nbs, seed = 72)

nfit = nfix
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
incs <- plotFit(fitW, df, fixedPar, type="SIgR", title = "GCT (n=12) -> GCT (n=12)")

fitPar = c(β = 0.2, D = 10, kappa = 2/9)
parEs <- onetimeParE(fitW, fitPar)

###################################################################### Different n

nfit = 6
fixedPar <- list(n = nfit, μ = μ, N = N, I0 = I0, ts = ts, nbs = nbs, T = T) 
fitW <- simplFit(startPar, fixedPar, df, sir.nll.g)
incd <- plotFit(fitW, df, fixedPar, type="SIgR", title = "GCT (n=12) -> GCT (n=6)")

parEd <- onetimeParE(fitW, fitPar)

row1 <- (incs | (parEs[[1]] | parEs[[2]] | parEs[[3]]))
row2 <- (incd | (parEd[[1]] | parEd[[2]] | parEd[[3]]))

print(row1 / row2)
