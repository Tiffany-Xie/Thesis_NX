library(deSolve)
library(devtools)
#install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(ggplot2) 
library(bbmle)
theme_set(theme_minimal())

######################################################################

simObs <- function(sinr, arp, nbs, seed) {
  set.seed(seed)
  inc <- diff(sinr[,"inc"]) # /ts
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  return(obs)
}

sir.nll.g <- function(obs, βe, De, n, μ, S0, I0, ts, T) {
  sir.null <- function(βe, De, obs) {
    out <- as.data.frame(SInRFlow(β=exp(βe), D=exp(De), n=n, μ=μ, S0=S0, I0=I0, ts=ts, T=T))
    nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=nbs, log=TRUE))
  }
}




###################################################################### Always work
β <- 0.4
D <- 10
#kappa <- 1/4
#fixn <- 12
n <- 4
μ <- 0.01

# Initial Value
S0 <- 999
I0 <- 1

# Time
ts <- 1
T <- 100

#sinner <- sinnerFlow(β, D, kappa, fixn, μ, ts, T)
sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
head(rowSums(sinr[,c(-1, -dim(sinr)[2])]))


arp <- 0.9
nbs <- 1000
inc <- diff(sinr[,"inc"]) /ts
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

sir.nll <- function(βe, De, obs){
  ## print(c(β=exp(βe), D=exp(De), n=4, μ=μ, S0, I0, ts, T))
  out <- as.data.frame(SInRFlow(β=exp(βe), D=exp(De), n=4, μ=μ, S0, I0, ts, T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=nbs, log=TRUE))
}

params0 <-list(βe=-0.5, De=2)
set.seed(33) # 33 work -0.4/0.6 X
fit0 <- mle2(sir.nll, start=params0, data=list(obs=obs))
print(fit0)
fit <- mle2(sir.nll, start=as.list(coef(fit0)), data=list(obs=obs))
print(fit)
p<-profile(fit)

plot(p, absVal=TRUE)
plot(timeSeq(ts, T), inc, type="l", xlab='Time, (Days)', ylab='I(t)')
plot(timeSeq(ts, T), obs, type="l", xlab='Time, (Days)', ylab='I(t)')

t <- timeSeq(ts, T)
mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(coef(fit)[["βe"]]), D=exp(coef(fit)[["De"]]), n=4, μ=0.01, S0, I0, ts, T)))
lines(diff(mod.prep$inc)~t, col = "red", lwd=3)

###################################################################### Change time
T <- 200

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
inc <- diff(sinr[,"inc"]) /ts
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

params0 <-list(βe=-0.5, De=2)
set.seed(33) # 33 work -0.4/0.6 X
fit0 <- invisible(mle2(sir.nll, start=params0, data=list(obs=obs))); fit0

###################################################################### Change N
T <- 100
S0 <- 99999

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
inc <- diff(sinr[,"inc"]) /ts
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

params0 <-list(βe=-0.5, De=2)
set.seed(33) # 33 work -0.4/0.6 X
fit0 <- invisible(mle2(sir.nll, start=params0, data=list(obs=obs))); fit0

###################################################################### Change time span
S0 <- 999
ts <- 0.1


sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
inc <- diff(sinr[,"inc"]) /ts
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

params0 <-list(βe=-0.5, De=2)
set.seed(33) # 33 work -0.4/0.6 X
fit0 <- invisible(mle2(sir.nll, start=params0, data=list(obs=obs))); fit0


###################################################################### Seeds - mle
seeds <- seq(1,100) 
β <- 0.4
D <- 10
n <- 4
μ <- 0.01
S0 <- 999
I0 <- 1
ts <- 1
T <- 100


succ_mle <- c()
sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
inc <- diff(sinr[,"inc"]) /ts
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

for (s in seeds) {
  set.seed(s)
  print("seed", s)
  params0 <-list(βe=-0.5, De=2)
  
  ans <- tryCatch({
    mle2(sir.nll, start=params0, data=list(obs=obs))
    1
  }, warning = function(w) {
    0.5
  }, error = function(e) {
    0
  }, finally = {
    1
  })
  
  succ_mle <- c(succ_mle, ans)
}

succ_mle

###################################################################### Seeds - NB
seeds <- seq(1,100) 
succ_nb <- c()
sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
inc <- diff(sinr[,"inc"]) /ts

for (s in seeds) {
  set.seed(s)
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  params0 <-list(βe=-0.5, De=2)
  
  ans <- tryCatch({
    mle2(sir.nll, start=params0, data=list(obs=obs))
    1
  }, warning = function(w) {
    0.5
  }, error = function(e) {
    0
  }, finally = {
    1
  })
  
  succ_nb <- c(succ_nb, ans)
}

succ_nb

###################################################################### Various initial value
β <- 0.4
D <- 10
n <- 4
μ <- 0.01
S0 <- 999
I0 <- 1
ts <- 1
T <- 100

arp <- 0.9
nbs <- 1000

βE=seq(-1, 1, 0.01); De=2


succ <- c()
sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
inc <- diff(sinr[,"inc"])
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

for (βe in βE) { 
  
  params0 <-list(βe=-βe, De=De)
  
  ans <- tryCatch({
    mle2(sir.nll, start=params0, data=list(obs=obs))
    1
  }, warning = function(w) {
    0.5
  }, error = function(e) {
    0
  }, finally = {
    1
  })
  
  succ <- c(succ, ans)
}

βE[which(succ==1)]

###################################################################### Large value problem
β <- 1
D <- 10
n <- 4
μ <- 0.01
S0 <- 999
I0 <- 1
ts <- 1
T <- 50

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
inc <- diff(sinr[,"inc"]) /ts
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

params0 <-list(βe=2, De=2)
set.seed(163) # 163 large value, but no error
fit0 <- invisible(mle2(sir.nll, start=params0, data=list(obs=obs))); fit0

# Even after using set.seed, it still sometimes generates error messages!!



quit()
######################################################################
# Try

tryCatch({
  
  log(-3)
  0
}, warning = function(w) {
  "A warning occurred"
  
}, error = function(e) {
  
  "An error occurred"
}, finally = {
  "This always runs"
})







