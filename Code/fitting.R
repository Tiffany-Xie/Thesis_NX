library(deSolve)
library(pseudoErlang)
library(gridExtra)
library(ggplot2) 
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
    
    dS <- μ*1000 - β*S*sum(I) - μ*S
    dI <- inflow - outflow
    dR <- inrate[n]*I[[n]] - μ*R
    inc <- β*S*sum(I)
    return(list(c(dS, dI, dR, inc)))
  })
}

######################################################################

SInRFlow <- function(β, mu, n, μ, ts, T) {
  gamma <- 1/mu
  outrate <- rep(n*gamma + μ, times=n)
  
  params <- list(β=β, n=n, μ=μ, outrate=outrate)
  states <- c(900, 100, numeric(n-1), 0, 0)
  names(states) <- c("S", paste0("I", 1:n), "R", "inc")
  return(Integration(params, states, ts, T, SIR))
}

sinnerFlow <- function(β, mu, kappa, n, μ, ts, T) {
  r <- kappa2r(kappa, n)
  a <- (1-1/r^n)/(mu*(1-1/r))
  #print(c(r, a))
  outrate <- a*r^(0:(n-1)) + μ
  
  params <- list(β=β, n=n, μ=μ, outrate=outrate)
  states <- c(9, 1, numeric(n-1), 0, 0)
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
ts <- 0.1
T <- 20

sinner <- sinnerFlow(β, mu, kappa, fixn, μ, ts, T)
sinr <- SInRFlow(β, mu, n, μ, ts, T)
head(rowSums(sinr[,c(-1, -dim(sinr)[2])]))

arp <- 0.9
inc <- diff(sinr[,"inc"])
nbs <- 100
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))


plot(timeSeq(ts, T), inc, type="l")

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
n=4
plot(obs~time,data=epi1, type='b', xlab='Day', ylab='I(t)')

sir.nll <- function(β, mu, μ){
  out <- as.data.frame(SInRFlow(β=β, mu=mu, n=4, μ=μ, ts=1, T=38))
  nll <- -sum(dpois(x=epi1$obs, lambda=tail(out$I4,15), log=TRUE))
}

params0 <-list(β=0.005, mu=7, μ=0.04)
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


