library(deSolve)
library(pseudoErlang)
library(gridExtra)
library(ggplot2) 
library(bbmle)
theme_set(theme_minimal())

######################################################################

Integration <- function(params, states, ts, T, model) {
  time <- timeSeq(ts, T, FALSE)
  soln <- ode(y = states,
              times = time,
              func = model,
              parms = params)
  return(soln)
}

######################################################################

SIR <- function(time, states, params) {
  with(params, {
    stopifnot(length(outrate)==n)
    S <- states[[1]]
    I <- states[2:(n+1)]
    R <- states[[n+2]]
    
    outflow <- outrate*I
    inrate <- outrate-μ
    inflow <- c(β*S*sum(I), (I*inrate)[1:(n-1)])
    
    dS <- μ*N - β*S*sum(I) - μ*S
    dI <- inflow - outflow
    dR <- inrate[n]*I[[n]] - μ*R
    inc <- β*S*sum(I)
    return(list(c(dS, dI, dR, inc)))
  })
}

######################################################################

SInRFlow <- function(β, mu, n, μ, S0, I0, ts, T) {
  gamma <- 1/mu
  outrate <- rep(n*gamma + μ, times=n)
  
  params <- list(β=β, n=n, μ=μ, N = S0+I0, outrate=outrate)
  states <- c(S0, I0, numeric(n-1), 0, 0)
  names(states) <- c("S", paste0("I", 1:n), "R", "inc")
  return(Integration(params, states, ts, T, SIR))
}

sinnerFlow <- function(β, mu, kappa, n, μ, S0, I0, ts, T) {
  r <- kappa2r(kappa, n)
  a <- (1-1/r^n)/(mu*(1-1/r))
  #print(c(r, a))
  outrate <- a*r^(0:(n-1)) + μ
  
  params <- list(β=β, n=n, μ=μ, N=S0+I0, outrate=outrate)
  states <- c(S0, I0, numeric(n-1), 0, 0)
  names(states) <- c("S", paste0("I", 1:n), "R", "inc")
  return(Integration(params, states, ts, T, SIR))
}

######################################################################

β <- 0.001
mu <- 10
kappa <- 1/4
fixn <- 12
n <- 4
μ <- 0.01
S0 <- 999
I0 <- 1
ts <- 0.1
T <- 30

sinner <- sinnerFlow(β, mu, kappa, fixn, μ, ts, T)
sinr <- SInRFlow(β, mu, n, μ, S0, I0, ts, T)
head(rowSums(sinr[,c(-1, -dim(sinr)[2])]))

arp <- 0.9
nbs <- 1000
inc <- diff(sinr[,"inc"])
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

plot(timeSeq(ts, T), inc, type="l", xlab='Time, (Days)', ylab='I(t)')
plot(timeSeq(ts, T), obs, type="l", xlab='Time, (Days)', ylab='I(t)')

#sir.nll <- function(βe, mue, μe, n=4){
#  S0 <- 999
#   I0 <- 1
#  out <- as.data.frame(SInRFlow(β=exp(βe), mu=exp(mue), n=n, μ=exp(μe), S0=S0, I0=I0, ts,T))
#  nll <- -sum(dpois(x=obs, lambda=diff(out$inc), log=TRUE))
#}

sir.nll <- function(β, mu, μ){
  S0 <- 999
  I0 <- 1
  out <- as.data.frame(SInRFlow(β=β, mu=mu, n=4, μ=μ, S0=S0, I0=I0, ts,T))
  nll <- -sum(dpois(x=obs, lambda=diff(out$inc), log=TRUE))
}

params0 <-list(β=0.04, mu=7, μ=0.006)
fit0 <- mle2(sir.nll, start=params0); fit0
fit <- mle2(sir.nll, start=as.list(coef(fit0))); fit
p<-profile(fit)


t <- timeSeq(ts, T)
mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=0.001346902, mu=15.000000137, n=4, μ=0.025714185, S0=999, I0=1, ts, T)))
lines(diff(mod.prep$inc)~t, col = "red")


# inc correct?
stop("Draft")
######################################################################
#library(remotes)
library(fitR)
#install_github("sbfnk/fitR", dependencies = TRUE)
library(bbmle)
data(epi)
ts = 0.1
T = 3
n=6
plot(obs~time,data=epi1, type='b', xlab='Day', ylab='I(t)')

sir.nll <- function(β, mu, μ){
  out <- as.data.frame(SInRFlow(β=β, mu=mu, n=6, μ=μ, ts=1, T=38))
  nll <- -sum(dpois(x=epi1$obs, lambda=out$I4, log=TRUE))
}

params0 <-list(β=0.005, mu=13, μ=0.004)
fit0 <- mle2(sir.nll, start=params0); fit0
fit <- mle2(sir.nll, start=as.list(coef(fit0))); fit
p<-profile(fit)

plot(p, absVal=TRUE)
plot(obs~time,data=epi1, type='b', xlab='Day', ylab='I(t)', col='red', ylim=c(0, 150))
t <- timeSeq(ts=1, T=38, FALSE)
mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=0.001194045, mu=3.590678810, n=4, μ=0.089549146, ts=1, T=38)))
lines(mod.prep$I4~t)
ylim=c(0, 150)


######################################################################
set.seed(123) 
n <- 1000 
size <- 10 
prob <- 0.5 
data <- rnbinom(n, size, prob)
data_df <- data.frame(counts = data)

ggplot(data_df, aes(x = counts)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Histogram of Negative Binomial Distribution",
       x = "Number of Failures",
       y = "Frequency")


