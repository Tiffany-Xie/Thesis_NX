
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
## Continuous data should be fitted with derlang and dperlang
gammaDelay <- function(n, mean, kappa){
  v <- rgamma(n, shape=1/kappa, scale=mean*kappa)
  print(c(mean = mean, kappa = kappa,
          Dmean = mean(v), Dkappa = sd(v)^2/mean(v)^2))
  return(v)
}

## Discrete data should be fitted another way (make boxed likelihood functions)
gammaDelayDiscrete <- function(n, mean, kappa){
  v <- round(rgamma(n, shape=1/kappa, scale=mean*kappa))
  print(c(mean = mean, kappa = kappa,
          Dmean = mean(v), Dkappa = sd(v)^2/mean(v)^2))
  return(v)
}

## Boxed likelihood example (not working yet!!!)
blgamma <- function(n){
  ## Do we need a special case for n=0?
  ## should we check that it's an integer and/or round?
  return(pgamma(n+1/2, log=TRUE) - pgamma(n-1/2, log=TRUE))
}

######################################################################

derlang <- function(x, box, rate, log=FALSE) { # box = n; rate = lambda
  density <- rate^box * x^(box-1) * exp(-rate*x) / factorial(box-1)
  if (log) {
    density <- log(density)
  }
  return(density)
} 

derlang2 <- function(x, D, box=3, log=FALSE) { # box = n; rate = lambda
  rate <- box * 1/D
  density <- rate^box * x^(box-1) * exp(-rate*x) / factorial(box-1)
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
  -sum(derlang(interval, box = round(1/exp(logkappa)), rate = 1/exp(logmean)*round(1/exp(logkappa)), log=TRUE))
}

ge.nll2 <- function(logmean, interval) {
  -sum(derlang2(interval, D = exp(logmean), log=TRUE))
}

gPE.nll <- function(logmean, logkappa, interval) {
  -sum(dperlang(interval, mean = exp(logmean), kappa = exp(logkappa), log=TRUE))
}

###################################################################### Pseudo Erlang

dperlang <- function(x, mean, kappa, box=12, log=FALSE) { # box=n; mean=D
  r <- kappa2r(kappa, box)
  a <- (1-1/r^box)/(mean*(1-1/r))
  density <- 0
  
  for (i in 1:box) {
    innerp <- 1
    
    for (j in 1:box) {
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

rperlang <- function(n, mean, kappa, box=12) {
  c <- 10 # gamma * 10  <<need to be improved>>
  samples <- numeric(n)
  i <- 1
  
  while (i<=n) {
    x <- rgamma(n=1, shape=1/kappa, scale=mean*kappa)
    u <- runif(1) # uniform(0,1)
    
    if (u <= dperlang(x, mean, kappa, box=box) / (c*dgamma(x, 1/kappa, scale=kappa*mean))) {
      samples[i] <- x
      i <- i + 1
    }
  }
  return(samples)
}

rperlang2 <- function(n, mean, kappa, box=12) {
  u <- runif(n)
  return(qperlang(u, mean, kappa, box=box))
}

pperlang <- function(x, mean, kappa, box=12, log=FALSE, offset=0) {
  r <- kappa2r(kappa, box)
  a <- (1-1/r^box)/(mean*(1-1/r))
  cdensity <- 0
  
  for (i in 1:box) {
    innerp <- 1
    
    for (j in 1:box) {
      if (j != i) {
        innerp <- innerp * r^(j-1) / (r^(j-1) - r^(i-1))
      }
    }
    cdensity <- cdensity + innerp * (1-exp(-x*a*r^(i-1)))
  }
  
  if (log) {
    cdensity <- log(cdensity-offset)
  }
  return(cdensity-offset)
}

## this would be much better if tseq were multiplied by mean (or mean*(1+kappa)) before pperlang
searchBound <- function(targetfx, mean, kappa, span=100, step=1, max=1000, box=12, log=FALSE) {
  unit <- mean*(1+kappa)
  span<-span*unit; step<-step*unit; max<-max*unit
  mintry <- 0
  while (mintry < max) {
    tseq <- seq(mintry, mintry+span, by=step)
    cdseq <- pperlang(tseq, mean, kappa)
    if (cdseq[length(cdseq)] >= targetfx) {
      break
    }
    mintry <- mintry+span
    if (mintry >= max) {
      return("You've reached the maximum")
    }
  }
  temp <- which(cdseq >= targetfx)[1]
  minindex <- temp - 1
  maxindex <- temp
  
  return(c(tseq[minindex]*0.9, tseq[maxindex]/0.9))
}

qperlang <- function(fx, mean, kappa, box=12, log=FALSE) {
  #if (fx < 1-1e-16)
  r <- kappa2r(kappa, box)
  a <- (1-1/r^box)/(mean*(1-1/r))
  
  f <- function(fxi) {
    bound <- searchBound(fxi, mean, kappa, box=box,log=log)
    #print(c(fxi, bound))
    root <- uniroot(pperlang, interval=c(bound[1], bound[2]), mean=mean, kappa=kappa, box=box, offset=fxi)
    root$root
  }
  
  x <- sapply(fx, f)
  return(x)
}

qperlang_o <- function(fx, mean, kappa, box=12, log=FALSE) {
  r <- kappa2r(kappa, box)
  a <- (1-1/r^box)/(mean*(1-1/r))
  
  f <- function(fxi) {
    upbound <- qgamma(fxi, box, rate=a)
    lowbound <- qgamma(fxi, box, rate=a*r^(box-1))
    root <- uniroot(pperlang, interval=c(lowbound, upbound), mean=mean, kappa=kappa, box=box, offset=fxi)
    root$root
  }
  
  x <- sapply(fx, f)
  return(x)
}
       
## need to find better xmax if possible (by formula / maybe it's not necessary..)
## R will round 1-1e-7 < # < 1 to 1
## <1-1e-15 qperlang no longer work; largest x = 120.8907
## qperlang can only take single fx, not vector  √fixed
## pperlang x > 60 gives 1

######################################################################

saveEnvironment()
