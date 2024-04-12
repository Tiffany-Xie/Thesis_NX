library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)
library(dplyr)

######################################################################
## data visualization
######################################################################

url <- "data/WHO-COVID-19-global-data.csv"
data <- read.csv(url)

data_c <- data %>%
  filter(Country == "China")

plot(seq(0,220), data_c$New_cases, type='l',xlab="Time, weeks", ylab="Incidence",
     main="Until Now (CHINA)")
plot(seq(0,220)[1:50], data_c$New_cases[1:50], type='l', xlab="Time, weeks", ylab="Incidence",
     main="First 50 Weeks (CHINA)")

######################################################################
## data visualization
######################################################################

nfit <- 20
startPar <- list(logβ=-1, logD=2, logkappa=-0.1)
fixedPar <- list(n = nfit, μ = 0.01, N = 12448000, I0 = 1, ts = 1, nbs = 1000, T = 30) 

sir.nll.g <- function(logβ, logD, logkappa, n, μ, N, I0, ts, T, nbs, obs) {
  
  out <- as.data.frame(sinnerFlow(β=exp(logβ), D=exp(logD), kappa=exp(logkappa), n=n, μ=μ, N=N, I0=I0, ts=ts, T=T))
  nll <- -sum(dnbinom(x=obs, mu=diff(out$Cinc), size=nbs, log=TRUE))
}

fit0 <- mle2(sir.nll.g, 
             data=list(obs=data_c$New_cases[1:30]),
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

df <- data.frame(Time = timeSeq(1,30), trueI = data_c$New_cases[1:30], fitI = diff(mod.prep$Cinc))

ggplot(df, aes(x=Time)) +
  geom_line(aes(y=trueI, color = "trueI"), linewidth=1) +
  geom_line(aes(y=fitI, color = "fitI")) +
  labs(x="Time, weeks", y = "Incidence")


