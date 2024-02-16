library(fitode)
library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)

#######################################################################

simObs <- function(sinr, arp, nbs, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  inc <- diff(sinr[,"Cinc"]) #/ts
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  df <- data.frame(times = sinr[,"time"],
                   #inc = inc,
                   obs = c(NA, obs))
  return(df)
}

#######################################################################

# generate SInR formula (list)
sinrFG <- function(n) {
  Ivec <- paste0("I",1:n)
  incidence <- sprintf("beta*S*(%s)/N", paste(Ivec, collapse = "+"))
  Itrans <- function(i, death=FALSE) {
    if (death) return(sprintf("(n*gamma + μ)*I%d", i))
    return(sprintf("n*gamma*I%d", i))
  }
  vars <- c("S", Ivec, "R", "Cinc")
  resp <- c(sprintf("μ*N - %s - μ*S", incidence),
            sprintf("%s - (n*gamma + μ)*I1", incidence),
            sprintf("%s - %s", Itrans((2:n)-1), Itrans(2:n, TRUE)),
            sprintf("%s - μ*R", Itrans(n)),
            incidence)
  
  fs <- purrr::map2(resp, vars, ~reformulate(.x, response = .y))
  fs <- lapply(fs, function(f) { environment(f) <- NULL; f })
  
  return(fs)
}

######################################################################

# generate initial value list
iniVG <- function(n) {
  Ivec <- paste0("I",1:n)
  vars <- c("S", Ivec, "R", "Cinc")
  inip <- c("N - i0",
            "i0",
            rep(0, n+1))
  
  Ini <- purrr::map2(inip, vars, ~reformulate(.x, response = .y))
  Ini <- lapply(Ini, function(f) { environment(f) <- NULL; f })
  
  return(Ini)
}

######################################################################

n <- 4

# include into a function ##
SIR_model
SIR_model <- odemodel(
  name="SInR (nbinom)",
  model=sinrFG(n),
  observation=list(
    obs ~ dnbinom(mu=Cinc, size=nbs)
  ),
  initial=iniVG(n),
  #diffnames="Cinc",
  par=c("beta", "gamma", "μ", "N", "i0", "nbs", "n"),
  #link=c(i0="logit")
)

######################################################################

β <- 0.2
D <- 10
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200

nbs <- 1000
arp <- 0.9

######################################################################

SIR_start <- c(beta=β, gamma=1/D, μ = μ, N=10000, i0=I0, nbs=nbs, n = n)
ss_SIR <- simulate(SIR_model,
                   parms=SIR_start, times=timeSeq(ts,T, FALSE))

######################################################################

sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)

df <- data.frame(Time = timeSeq(ts, T),
                 currentInc = diff(sinr[,'Cinc']),
                 fitodeInc = diff(ss_SIR[,'Cinc']))

######################################################################

print(ggplot(df, aes(x=Time))
      + geom_line(aes(y=currentInc, color="original inc"), linewidth = 1.5)
      + geom_line(aes(y=fitodeInc, color="fitode inc"), linetype = "dashed", linewidth = 2) 
      + labs(y = "Incidence", title = "Original SInR Inc vs. FitOde SInR Inc (n=4)")
)

###################################################################### Fit

n <- 4

# include into a function ##
SIR_model <- odemodel(
  name="SInR (nbinom)",
  model=sinrFG(n),
  observation=list(
    obs ~ dnbinom(mu=Cinc, size=nbs)
  ),
  initial=iniVG(n),
  diffnames="Cinc",
  par=c("beta", "gamma", "μ", "N", "i0", "nbs", "n"),
  #link=c(i0="logit")
)

SIR_start <- c(beta=β, gamma=1/D, μ = μ, N=10000, i0=I0, nbs=nbs, n = n)
ss_SIR <- simulate(SIR_model,
                   parms=SIR_start, times=timeSeq(ts,T, FALSE))

######################################################################
# Fit
df = simObs(sinr, arp, nbs, seed = 73)

SIR_fit <- fitode(
  SIR_model,
  data=df$obs,
  start=SIR_start
)

####### ######## ######### 
SierraLeone2014b <- rbind(
  c(times=SierraLeone2014$times[1] -
      diff(SierraLeone2014$times)[1], confirmed=NA),
  SierraLeone2014
)

SIR_model <- odemodel(
  name="SIR (nbinom)",
  model=list(
    S ~ - beta * S * I/N,
    I ~ beta * S * I/N - gamma * I,
    R ~ gamma * I
  ),
  observation=list(
    confirmed ~ dnbinom(mu=R, size=phi)
  ),
  initial=list(
    S ~ N * (1 - i0),
    I ~ N * i0,
    R ~ 0
  ),
  diffnames="R",
  par=c("beta", "gamma", "N", "i0", "phi"),
  link=c(i0="logit")
)

SIR_start <- c(beta=70, gamma=60, N=40000, i0=0.0004, phi=6)
ss_SIR <- simulate(SIR_model,
                   parms=SIR_start, times=SierraLeone2014b$times)
plot(SierraLeone2014)
lines(ss_SIR$times, ss_SIR$R)

SIR_fit <- fitode(
  SIR_model,
  data=SierraLeone2014b,
  start=SIR_start
)

summary(SIR_fit)



