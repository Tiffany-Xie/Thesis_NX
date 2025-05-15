library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 12))
library(bbmle)
library(pseudoErlang)
library(patchwork)

library(shellpipes)
startGraphics(height=3, width=11)

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
kappa<- 0.3

######################################################################

parEst <- function(nfix, cnfit, fitPar, seeds) {
  sigr <- sinnerFlow(β, D, kappa, nfix, μ, N, I0, ts, T)
  startPar_g <- list(logβ=-1, logD=2, logkappa=-1)
  startPar_l <- list(logβ=-1, logD=2)
  
  plots <- list()
  for (i in names(fitPar)) {
    data_g <- data.frame(variable = integer(),
                         estimate = numeric(),
                         lower = numeric(),
                         upper = numeric())
    data_l <- data.frame(variable = integer(),
                         estimate = numeric(),
                         lower = numeric(),
                         upper = numeric())
    for (n in cnfit) {
      for (seed in seeds) {
        df <- simObs(sigr, arp, nbs, seed = seed)
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
        else if (i == "kappa") {
          data_l <- rbind(data_l, data.frame(variable=n, estimate=1/n, lower=1/n, upper=1/n)) 
        }
      }
    }
    print(data_l)
    #print(data_g)
    #print(ggplot(data_g, aes(x=variable, y=estimate)) + 
    #  geom_hline(yintercept=fitPar[[i]], col='red', linewidth=1) +
    #  geom_point() +
    #  coord_flip()+
    #  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
    #  geom_point(data=data_l, aes(x=variable, y=estimate),
    #             shape=23, color="black", fill="yellow2", size=3) +
    #  geom_errorbar(data=data_l, aes(ymin=lower, ymax=upper), width=0.2, linewidth=1.5, color="gray70", alpha=0.8) +
    #  #scale_shape_manual(values = c(20, 5)) + 
    #  #scale_color_manual(values = c("black", "yellow2"), labels=c("GCT", "LCT")) +
    #  labs(x="Substage Numbers", y = paste0("Estimate ", ifelse(i == "β", "beta", i)))
    #  )
    #plots <- c(plots, list(p))
    data_l$Model <- 'data_l'
    data_g$Model <- 'data_g'
    combined_data <- rbind(data_l, data_g)
    p <- ggplot(combined_data, aes(x = variable, y = estimate, color=Model)) +
            geom_point(size=2) +
            geom_hline(yintercept=fitPar[[i]], col='red', linewidth=1) +
            geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
            #scale_shape_manual(values = c(20, 5)) +
            scale_color_manual(values = c("black", "green2"), labels=c("GCT", "LCT")) +
            labs(x = "Substage Numbers", y = paste0("Estimate ", ifelse(i == "β", "beta", i))) +
            coord_flip()
    plots <- c(plots, list(p))
  }
  return(plots)
}

nfix <- 4
cnfit <- seq(4,12, by=2)
fitPar = c(β=0.2, D=10, kappa=0.3)
seeds <- c(71, 77, 121, 123)

plots <- parEst(nfix, cnfit, fitPar, seeds)

######################################################################

plots[[1]] | plots[[2]] | plots[[3]]


