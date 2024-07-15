
library(shellpipes)

simObs <- function(sinr, arp, nbs, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  inc <- diff(sinr[,"Cinc"]) #/ts
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  df <- data.frame(Time = (sinr[,"time"][-1] + sinr[,"time"][-dim(sinr)[1]])/2,
                   inc = inc,
                   obs = obs)
  return(df)
}

######################################################################

## likelihood function with parameter-trajectory tracing
sir.nll <- function(logβ, logD, n, μ, N, I0, ts, T, nbs, obs) {
  trace_logbeta <<- c(trace_logbeta, logβ)
  trace_logD <<- c(trace_logD, logD)
  
  out <- as.data.frame(SInRFlow(β=exp(logβ), D=exp(logD), n=n, μ=μ, N=N, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

## Similar function for geometric pseudo-Erlang (different params)
sir.nll.g <- function(logβ, logD, logkappa, n, μ, N, I0, ts, T, nbs, obs) {
  trace_logbeta <<- c(trace_logbeta, logβ)
  trace_logD <<- c(trace_logD, logD)
  trace_logkappa <<- c(trace_logkappa, logkappa)
  
  out <- as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa), n=n, μ=μ, N=N, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

######################################################################

simplFit <- function(startPar, fixedPar, datadf, likelihood.m, optMethod="Nelder-Mead") {
  obs <- round((datadf$obs)/arp)
  trace_logbeta <<- c()
  trace_logD <<- c()
  trace_logkappa <<- c()
  
  fit0 <- mle2(likelihood.m, 
               data=list(obs=obs),
               start=startPar, 
               fixed=fixedPar,
               method=optMethod,
               control=list(trace = 0))
  
  return(list(fit=fit0, trace_logbeta=trace_logbeta, trace_logD=trace_logD))
}

######################################################################

plotFit <- function(fitW, df, fPar, type="SInR", title = "Fitting Result") {
  logβ <- coef(fitW$fit)[["logβ"]]
  logD <- coef(fitW$fit)[["logD"]]
  
  if (type=="SInR") {
    mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(logβ), D=exp(logD), 
                                                     n=fPar$n, μ=fPar$μ, N=fPar$N, 
                                                     I0=fPar$I0, ts=fPar$ts, T=fPar$T)))
  }
  else if (type=="SIgR") {
    logkappa <- coef(fitW$fit)[["logkappa"]]
    mod.prep <- as.data.frame(as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa),
                                                       n=fPar$n, μ=fPar$μ, N=fPar$N, 
                                                       I0=fPar$I0, ts=fPar$ts, T=fPar$T)))
  }
  
  df["fitInc"] = diff(mod.prep$Cinc)
  #mse <- mean((df$inc - df$fitInc)^2)
  RMSLE <- sqrt(sum(log((df$fitInc+1)/(df$inc+1))^2)/dim(df)[1])
  
  p <- ggplot(df, aes(x=Time)) +
    geom_line(aes(y=obs, color='Observed')) +
    geom_line(aes(y=inc, color = 'Incidence'), linewidth=1) +
    geom_line(aes(y=fitInc, color = 'Fit Incidence'), linewidth=1.5, alpha=0.8) +
    labs(title = paste0(title, " | RMSLE = ", round(RMSLE*100, digits=2), "%"), x = "Time, (days)", y = "Incidence")
  return(p)
}


plotTrace <- function(fitW) {
  trace_logbeta <- fitW$trace_logbeta
  trace_logD <- fitW$trace_logD
  
  pardf <- data.frame(order = seq(0:(length(trace_logbeta)-1)), 
                      logbeta = trace_logbeta,
                      logD = trace_logD)
  
  ggplot(pardf, aes(x = logD, y = logbeta)) +
    geom_point(aes(color = order)) +
    geom_path(alpha = 0.5) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(x = "log(D)", y = "log(beta)", color = "Time")
}

######################################################################
######################################################################

derlang <- function(x, shape, rate, log=FALSE) { # shape = n; rate = lambda
  density <- rate^shape * x^(shape-1) * exp(-rate*x) / factorial(shape-1)
  if (log) {
    density <- log(density)
  }
  return(density)
} 

derlang2 <- function(x, D, shape=3, log=FALSE) { # shape = n; rate = lambda
  rate <- shape * 1/D
  density <- rate^shape * x^(shape-1) * exp(-rate*x) / factorial(shape-1)
  if (log) {
    density <- log(density)
  }
  return(density)
}

dperlang <- function(x, mean, kappa, shape=12, log=FALSE) { # shape=n; mean=D
  r <- kappa2r(kappa, shape)
  a <- (1-1/r^shape)/(mean*(1-1/r))
  density <- 0
  
  for (i in 1:shape) {
    innerp <- 1
    
    for (j in 1:shape) {
      if (j != i) {
        innerp <- innerp * r^(j-1)/(r^(j-1) - r^(i-1))
      }
    }
    density <- density + innerp * a * r^(i-1) * exp(-x*a*r^(i-1))
  }
  
  if (log) {
    density <- log(density)
  }
  return(density)
} 

######################################################################

gg.nll <- function(logmean, logkappa, interval) {
  -sum(dgamma(interval, shape = 1/exp(logkappa), scale = exp(logmean+logkappa), log=TRUE))
}

ge.nll <- function(logmean, logkappa, interval) {
  -sum(derlang(interval, shape = round(1/exp(logkappa)), rate = 1/exp(logmean)*round(1/exp(logkappa)), log=TRUE))
}

ge.nll2 <- function(logmean, interval) {
  -sum(derlang2(interval, D = exp(logmean), log=TRUE))
}

gPE.nll <- function(logmean, logkappa, interval) {
  -sum(dperlang(interval, mean = exp(logmean), kappa = exp(logkappa), log=TRUE))
}

######################################################################

rperlang <- function(n, mean, kappa, shape=12) {
  c <- 10 # gamma * 10  <<need to be improved>>
  samples <- numeric(n)
  i <- 1
  
  while (i<=n) {
    x <- rgamma(n=1, shape=1/kappa, scale=mean*kappa)
    u <- runif(1) # uniform(0,1)
    
    if (u <= dperlang(x, mean, kappa, shape=shape) / (c*dgamma(x, 1/kappa, scale=kappa*mean))) {
      samples[i] <- x
      i <- i + 1
    }
  }
  return(samples)
}

######################################################################

saveEnvironment()

