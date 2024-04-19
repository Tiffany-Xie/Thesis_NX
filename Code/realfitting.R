library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)
library(dplyr)

######################################################################
## data visualization
######################################################################

path <- "https://raw.githubusercontent.com/Tiffany-Xie/Thesis_NX/main/data/WHO_covid.csv?token=GHSAT0AAAAAACN7IGL3P7DCW6BOB7ZBLLDAZRABGEQ"
data <- read.csv(path)

data_c <- data %>% 
  filter(Country == "China")

plot(seq(0,dim(data_c)[1]-1), data_c$New_cases, type='l', xlab="Time, weeks", ylab="Incidence",
     main="Until Now (CHINA)")
plot(seq(0,dim(data_c)[1]-1)[1:50], data_c$New_cases[1:50], type='l', xlab="Time, weeks", ylab="Incidence",
     main="First 50 Weeks (CHINA)")
plot(seq(0,dim(data_c)[1]-1)[150:175], data_c$New_cases[150:175], type='l', xlab="Time, weeks", ylab="Incidence",
     main="Maximum Peak (CHINA)")

######################################################################
## fitting
######################################################################

nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = 0.01, N = 1.4097e9, I0 = 271745, ts = 1, nbs = 1000, T = 26) # 171745

sir.nll.g <- function(logβ, logD, logkappa, n, μ, N, I0, ts, T, nbs, obs) {
  
  out <- as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa), n=n, μ=μ, N=N, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

fit0 <- mle2(sir.nll.g, 
             data=list(obs=data_c$New_cases[150:175]),
             start=startPar, 
             fixed=fixedPar,
             method="Nelder-Mead",
             control=list(trace = 0))

logβ <- coef(fit0)[["logβ"]]
logD <- coef(fit0)[["logD"]]
logkappa <- coef(fit0)[["logkappa"]]
mod.prep <- as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa),
                                     n=fixedPar$n, μ=fixedPar$μ, N=fixedPar$N, 
                                     I0=fixedPar$I0, ts=fixedPar$ts, T=fixedPar$T))

df <- data.frame(Time = timeSeq(1,26), trueI = data_c$New_cases[150:175], fitI = diff(mod.prep$Cinc))

ggplot(df, aes(x=Time)) +
  geom_line(aes(y=trueI, color = "trueI"), linewidth=1) +
  geom_line(aes(y=fitI, color = "fitI")) +
  labs(x="Time, weeks", y = "Incidence")

######################################################################
## Measles
######################################################################

measles_url <- "https://raw.githubusercontent.com/mac-theobio/Disease_data/master/outputs/ewmeas.ssv"
data <- read.table(measles_url, header=FALSE, sep=" ")
names(data) = c("Time", "Cases")

which.min(data$Cases[330:450])
plot(data$Time, data$Cases, type="l", xlab="Time", ylab="New Cases")
data <- as.data.frame(data)
ggplot(data, aes(x=Time, y=Cases)) + 
  geom_line() +
  geom_point(data=data.frame(x=data$Time[c(349, 454)], y=data$Cases[c(349, 454)]), aes(x,y), color="red",size=3)

which.min(subset(data, Time>1956.5 & Time<1957.5)$Cases)
subset(data, Time>1956.5 & Time<1957.5)[11,]
onepeak <- data[349:454,] #299:402
ggplot(onepeak, aes(x=Time, y=Cases)) + geom_line()


nfit <- 12
startPar <- list(logβ=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, μ = 0.01, N = 160180000, I0 = 1, ts = 7, nbs = 1000, T = 742) # 171745

sir.nll.g <- function(logβ, logD, logkappa, n, μ, N, I0, ts, T, nbs, obs) {
  
  out <- as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa), n=n, μ=μ, N=N, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

fit0 <- mle2(sir.nll.g, 
             data=list(obs=onepeak$Cases),
             start=startPar, 
             fixed=fixedPar,
             method="Nelder-Mead",
             control=list(trace = 0))

logβ <- coef(fit0)[["logβ"]]
logD <- coef(fit0)[["logD"]]
logkappa <- coef(fit0)[["logkappa"]]
mod.prep <- as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa),
                                     n=fixedPar$n, μ=fixedPar$μ, N=fixedPar$N, 
                                     I0=fixedPar$I0, ts=fixedPar$ts, T=fixedPar$T))

df <- data.frame(Time = timeSeq(7/365,7*106/365)+1954.696, trueI = onepeak$Cases, fitI = diff(mod.prep$Cinc))

ggplot(df, aes(x=Time)) +
  geom_line(aes(y=trueI, color = "trueI"), linewidth=1) +
  geom_line(aes(y=fitI, color = "fitI")) +
  labs(x="Time, weeks", y = "Incidence")

######################################################################
## Measles seasonal 
######################################################################

SIR <- function(time, states, params) {
  with(params, {
    stopifnot(length(outrate)==n)
    β <- β0*(1 + sigma*sinpi(2/P*time))
    
    S <- states[[1]]
    I <- states[2:(n+1)]
    R <- states[[n+2]]
    
    Cinc <- β*S*sum(I)/N
    
    outflow <- outrate*I
    inrate <- outrate-μ
    inflow <- c(Cinc, (I*inrate)[0:(n-1)])
    
    dS <- μ*N - Cinc - μ*S
    dI <- inflow - outflow
    dR <- inrate[n]*I[[n]] - μ*R
    dC <- Cinc
    
    return(list(c(dS, dI, dR, dC)))
  })
}

sinnerFlow <- function(β0, sigma, P, D, kappa, n, μ, N, I0, ts, T) {
  if (kappa == 1/n) {
    r = n/D
    a = 1
  }
  else {
    r <- kappa2r(kappa, n)
    a <- (1-1/r^n)/(D*(1-1/r))
  }
  #print(c(r, a))
  outrate <- a*r^(0:(n-1)) + μ
  
  params <- list(β0=β0, sigma=sigma, P=P, n=n, μ=μ, N=N, outrate=outrate)
  states <- c(N-I0, I0, numeric(n-1), 0, 0)
  names(states) <- c("S", paste0("I", 1:n), "R", "Cinc")
  return(Integration2(params, states, ts, T, SIR))
}
###
nfit <- 12
startPar <- list(logβ0=-1, logD=2, logkappa=-0.2)
fixedPar <- list(n = nfit, sigma=0.4, P=365*2, μ = 0.01, N = 160180000, I0 = 3700, ts = 7, nbs = 1000, T = 7*991) # 171745

sir.nll.g <- function(logβ0, logD, logkappa, n, sigma, P, μ, N, I0, ts, T, nbs, obs) {
  
  out <- as.data.frame(sinnerFlow(β=exp(logβ0), D=exp(logD), kappa=exp(logkappa), sigma=sigma, P=P, n=n, μ=μ, N=N, I0=I0, ts=ts, T=6937))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

fit0 <- mle2(sir.nll.g, 
             data=list(obs=data$Cases),
             start=startPar, 
             fixed=fixedPar,
             method="Nelder-Mead",
             control=list(trace = 0))

logβ0 <- coef(fit0)[["logβ0"]]
logD <- coef(fit0)[["logD"]]
logkappa <- coef(fit0)[["logkappa"]]
mod.prep <- as.data.frame(sinnerFlow(β0=exp(logβ0), D=exp(logD), kappa=exp(logkappa),
                                     n=fixedPar$n, μ=fixedPar$μ, N=fixedPar$N, 
                                     I0=fixedPar$I0, ts=fixedPar$ts, T=fixedPar$T, 
                                     sigma=fixedPar$sigma, P=fixedPar$P))

df <- data.frame(Time = timeSeq(7/365,7*991/365)+1948.027, trueI = data$Cases, fitI = diff(mod.prep$Cinc))

ggplot(df, aes(x=Time)) +
  geom_line(aes(y=trueI, color = "trueI"), linewidth=1) +
  geom_line(aes(y=fitI, color = "fitI")) +
  labs(x="Time, weeks", y = "Incidence", title = "Whole Data Fitting")


