library(deSolve)
library(pseudoErlang)
library(bbmle)

set.seed(43) ## to get an error
set.seed(42)

Integration <- function(params, states, ts, T, model) {
	print(params)
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
# Parameters
β <- 0.001
mu <- 10
#kappa <- 1/4
#fixn <- 12
n <- 4
μ <- 0.01

# Initial Value
S0 <- 999
I0 <- 1

# Time
ts <- 0.1
T <- 30

#sinner <- sinnerFlow(β, mu, kappa, fixn, μ, ts, T)
sinr <- SInRFlow(β, mu, n, μ, S0, I0, ts, T)
head(rowSums(sinr[,c(-1, -dim(sinr)[2])]))
head(sinr)

arp <- 0.9
nbs <- 1000
inc <- diff(sinr[,"inc"])
obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))

print(obs)

sir.nll <- function(βe, mue, obs){
  out <- as.data.frame(SInRFlow(β=exp(βe), mu=exp(mue), n=4, μ=μ, S0=999, I0=1, ts, T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=nbs, log=TRUE))
}

params0 <-list(βe=-6, mue=1)
fit0 <- mle2(sir.nll, start=params0, data=list(obs=obs)); fit0
fit <- mle2(sir.nll, start=as.list(coef(fit0)), data=list(obs=obs)); fit
p<-profile(fit)

plot(p, absVal=TRUE)
plot(timeSeq(ts, T), inc, type="l", xlab='Time, (Days)', ylab='I(t)')
plot(timeSeq(ts, T), obs, type="l", xlab='Time, (Days)', ylab='I(t)')

t <- timeSeq(ts, T)
mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(coef(fit)[["βe"]]), mu=exp(coef(fit)[["mue"]]), n=4, μ=0.01, S0=999, I0=1, ts, T)))
lines(diff(mod.prep$inc)~t, col = "red", lwd=3)

# change time
# change to exponential 

# inc correct?
# change time/ts influence result

quit()

######################################################################
load(system.file("vignetteData","orob1.rda",package="bbmle"))
summary(orob1)

X <- model.matrix(~dilution, data = orob1)
ML1 <- function(prob1,prob2,prob3,theta,x) {
  prob <- c(prob1,prob2,prob3)[as.numeric(x$dilution)]
  size <- x$n
  -sum(dbetabinom(x$m,prob,size,theta,log=TRUE))
}


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


