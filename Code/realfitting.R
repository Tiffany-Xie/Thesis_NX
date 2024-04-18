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
  geom_point(data=data.frame(x=data$Time[c(299,402)], y=data$Cases[c(299,402)]), aes(x,y), color="red",size=3)


