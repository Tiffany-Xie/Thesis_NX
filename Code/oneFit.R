library(deSolve)
library(devtools)
library(pseudoErlang)
library(ggplot2) 
library(bbmle)

theme_set(theme_minimal())

sir.nll <- function(βe, De, obs, n, μ, S0, I0, ts, T){
	## print(c(β=exp(βe), D=exp(De), n=4, μ=μ, S0, I0, ts, T))
  out <- as.data.frame(SInRFlow(β=exp(βe), D=exp(De), n=n, μ=μ, S0=S0, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=nbs, log=TRUE))
}

simObs <- function(β, D, n, μ, S0, I0, ts, T, rp, nbsize){
	sim <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
	inc <- diff(sim[,"inc"])/ts
	return(rnbinom(mu=rp*inc, size=nbs, n=length(inc)))
}

simpleFit <- function(startPar, fixedPar, inc, seed){
	set.seed(seed)
	fit0 <- mle2(sir.nll, start=startPar, fixed=fixedPar, data=list(obs=obs))
}

β <- 0.4
D <- 10
n <- 4
μ <- 0.01
S0 <- 999
I0 <- 1
ts <- 1
T <- 100
rp <- 0.9
nbs <- 1000
seed <- 33

startPar <-list(βe=-0.5, De=2)
fixedPar <- list(
	n = n, μ = μ, S0 = S0, I0 = I0, ts = ts, T = T
)

obs <- simObs(β, D, n, μ, S0, I0, ts, T, rp, nbsize)
summary(obs)

print(simpleFit(startPar, fixedPar, inc, seed))
