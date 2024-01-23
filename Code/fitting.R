library(deSolve)
library(pseudoErlang)
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
    
    inc <- β*S*sum(I)/N
    
    outflow <- outrate*I
    inrate <- outrate-μ
    inflow <- c(inc, (I*inrate)[1:(n-1)])
    
    dS <- μ*N - inc - μ*S
    dI <- inflow - outflow
    dR <- inrate[n]*I[[n]] - μ*R
    dC <- inc
    
    return(list(c(dS, dI, dR, dC)))
  })
}

######################################################################

SInRFlow <- function(β, D, n, μ, S0, I0, ts, T) {
  gamma <- 1/D
  outrate <- rep(n*gamma + μ, times=n)
  
  params <- list(β=β, n=n, μ=μ, N = S0+I0, outrate=outrate)
  states <- c(S0, I0, numeric(n-1), 0, 0)
  names(states) <- c("S", paste0("I", 1:n), "R", "inc")
  return(Integration(params, states, ts, T, SIR))
}

sinnerFlow <- function(β, D, kappa, n, μ, S0, I0, ts, T) {
  r <- kappa2r(kappa, n)
  a <- (1-1/r^n)/(D*(1-1/r))
  #print(c(r, a))
  outrate <- a*r^(0:(n-1)) + μ
  
  params <- list(β=β, n=n, μ=μ, N=S0+I0, outrate=outrate)
  states <- c(S0, I0, numeric(n-1), 0, 0)
  names(states) <- c("S", paste0("I", 1:n), "R", "inc")
  return(Integration(params, states, ts, T, SIR))
}

######################################################################
# Parameters
β <- 0.3
D <- 10
#kappa <- 1/4
#fixn <- 12
n <- 4
μ <- 0.01

# Initial Value
S0 <- 9999
I0 <- 1

# Time
ts <- 1
T <- 100

#sinner <- sinnerFlow(β, D, kappa, fixn, μ, ts, T)
sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
head(rowSums(sinr[,c(-1, -dim(sinr)[2])]))

arp <- 0.9
nbs <- 1000
inc <- diff(sinr[,"inc"])
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

sir.nll <- function(βe, De, obs){
  out <- as.data.frame(SInRFlow(β=exp(βe), D=exp(De), n=4, μ=μ, S0, I0, ts, T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=nbs, log=TRUE))
}

params0 <-list(βe=-1, De=1)
fit0 <- mle2(sir.nll, start=params0, data=list(obs=obs)); fit0
fit <- mle2(sir.nll, start=as.list(coef(fit0)), data=list(obs=obs)); fit
p<-profile(fit)

plot(p, absVal=TRUE)
plot(timeSeq(ts, T), inc, type="l", xlab='Time, (Days)', ylab='I(t)')
plot(timeSeq(ts, T), obs, type="l", xlab='Time, (Days)', ylab='I(t)')

t <- timeSeq(ts, T)
mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(coef(fit)[["βe"]]), D=exp(coef(fit)[["De"]]), n=4, μ=0.01, S0, I0, ts, T)))
lines(diff(mod.prep$inc)~t, col = "red", lwd=3)

quit()
######################################################################
seeds <- seq(1,50)
β <- 0.001
D <- 10
n <- 4
μ <- 0.01
S0 <- 999
I0 <- 1
ts <- 0.1
T <- 50

arp <- 0.9
nbs <- 1000

succ <- c()

for (s in seeds) {
  set.seed(s)
  sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
  
  inc <- diff(sinr[,"inc"])
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  
  params0 <-list(βe=-6, De=1)
  
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





