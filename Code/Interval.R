library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)

library(shellpipes)
loadEnvironments()
# source("./tempfunc.R")
######################################################################

## Continuous data should be fitted with derlang and dperlang
gammaDelay <- function(n, mean, kappa){
  v <- rgamma(n, shape=1/kappa, scale=mean*kappa)
  print(c(mean = mean, kappa = kappa,
          Dmean = mean(v), Dkappa = sd(v)^2/mean(v)^2))
  return(v)
}

## Discrete data should be fitted another way (make boxed likelihood functions)
gammaDelayDiscrete <- function(n, mean, kappa){
  v <- round(rgamma(n, shape=1/kappa, scale=mean*kappa))
  print(c(mean = mean, kappa = kappa,
          Dmean = mean(v), Dkappa = sd(v)^2/mean(v)^2))
  return(v)
}

## Boxed likelihood example (not working yet!!!)
blgamma <- function(n){
	## Do we need a special case for n=0?
	## should we check that it's an integer and/or round?
	return(pgamma(n+1/2, log=TRUE) - pgamma(n-1/2, log=TRUE))
}

######################################################################
## Gamma -> Erlang
######################################################################

## k = 1/shape
## D = shape*scale = scale/kappa

## start/end times
## start is a uniform
## end is start + v [above]
## round them both before subtracting

set.seed(223)
mean = 7
kappa = 0.3
n = 2000

g_interval = gammaDelay(n,mean,kappa) # v
time = seq(0, max(g_interval), length.out=n)

startPar = list(logmean = 2.5, logkappa = -2) # 2 -1
fit <- mle2(ge.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 fit_gamma = derlang(time, box=round(1/fitkappa), rate=1/fitmean*round(1/fitkappa)))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Erlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = paste0("gamma -> Erlang |", 
                                                   " fitmean=", round(fitmean, 3), 
                                                   ", fitkappa=", round(fitkappa, 3),
                                                   ", Loglik=", round(logLik(fit), 3)))

######################################################################
## Gamma -> Pseudo Erlang
######################################################################

startPar = list(logmean = 2, logkappa = -1)
fit <- mle2(gPE.nll, 
            data = list(interval = g_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = g_interval,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 fit_gamma = dperlang(time, mean=fitmean, kappa=fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Pseudo Erlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = paste0("gamma -> Pseudo Erlang |", 
                                                   " fitmean=", round(fitmean, 3), 
                                                   ", fitkappa=", round(fitkappa, 3),
                                                   ", Loglik=", round(logLik(fit), 3)))
# actual data kappa
p <- ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density)), bins=45)
hist_data <- ggplot_build(p)$data[[1]]

actmean = print(sum(g_interval)/length(g_interval))
span = (hist_data$xmax[1] + hist_data$xmin[1])/2
actvar = sum(hist_data$x^2 * hist_data$density * span) - actmean^2
actkappa = print(actvar/actmean^2)

######################################################################
## Inverse Pseudo Erlang CDF
######################################################################
mean=7
kappa=0.3
cd <- seq(0.001, 0.999, 0.001)
system.time(inversed_time_pe <- qperlang(cd, mean, kappa))
system.time(inversed_time_pe <- qperlang_o(cd, mean, kappa))
inversed_time_g <- qgamma(cd, 1/kappa, scale=mean*kappa)
df <- data.frame(CD = cd, 
                 Interval_pe = inversed_time_pe,
                 Interval_g = inversed_time_g)

ggplot(df, aes(x=CD)) +
  geom_line(aes(y=Interval_pe^2, color='Pseudo Erlang(sq)'), linewidth=1) +
  geom_line(aes(y=Interval_g^2, color='Gamma(sq)'), linewidth=1)


plot(x=cd, y=inversed_time, type="l",
     xlab="Cumulative Density",
     ylab="Interval",
     main="Inversed Pseudo Erlang CDF")

#####################################################################
## Random Pseudo Erlang numbers
######################################################################

set.seed(1001)
n <- 1000
mean <- 7
kappa <- 0.3

pe_interval <- rperlang2(n, mean, kappa) # <<slow if c->inf>>
time = seq(0, max(pe_interval), length.out=n)
df <- data.frame(Time=time, interval=pe_interval,
                 perlang=dperlang(time, mean, kappa),
                 gamma=dgamma(time, 1/kappa, scale=mean*kappa))
p <- ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density)))+
  geom_line(aes(x=Time, y=perlang, color="PErlang"), linewidth=1.5) +
  geom_line(aes(x=Time, y=gamma, color="Gamma"), linewidth=1.5) +
  labs(title="Random PErlang (inversion-based)")

hist_data <- ggplot_build(p)$data[[1]]
span <- hist_data$x[2] + hist_data$x[1]
print(sum(hist_data$density * 1.111401))

## manually
bins <- 30
binwidth <- (max(pe_interval) - min(pe_interval))/bins
breaks <- seq(min(pe_interval), max(pe_interval), by=binwidth)
binned <- cut(pe_interval, breaks=breaks, include.lowest = TRUE, right = TRUE)
counts <- table(binned)
counts_df <- as.data.frame(counts)
colnames(counts_df) <- c("Days", "Count")

counts_df$Density <- counts_df$Count/(length(pe_interval)*binwidth) 
counts_df$midpoint <- seq(min(pe_interval)+binwidth/2, max(pe_interval), by=binwidth)
ggplot(counts_df, aes(x=midpoint, y=Density)) +
         geom_bar(stat = "identity", width = 1.08)

sum(binwidth*counts_df$Density)
######################################################################
## Pseudo Erlang -> Pseudo Erlang
######################################################################

startPar = list(logmean = 2.5, logkappa = -0.5)

fit <- mle2(gPE.nll, 
            data = list(interval = pe_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = pe_interval,
                 perlang = dperlang(time, mean, kappa),
                 fit_perlang = dperlang(time, fitmean, fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=perlang, color="PErlang"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_perlang, color="PErlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = paste0("PErlang -> PErlang |", 
                                                   " fitmean=", round(fitmean, 3), 
                                                   ", fitkappa=", round(fitkappa, 3),
                                                   ", Loglik=", round(logLik(fit), 3)))

######################################################################
## Pseudo Erlang -> Gamma
######################################################################

startPar = list(logmean = 2.5, logkappa = -0.5)

fit <- mle2(gg.nll, 
            data = list(interval = pe_interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = pe_interval,
                 perlang = dperlang(time, mean, kappa),
                 fit_gamma = dgamma(time, 1/fitkappa, scale=fitmean*fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=perlang, color="PErlang"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_gamma, color="Gamma after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = paste0("PErlang -> gamma |", 
                                                   " fitmean=", round(fitmean, 3), 
                                                   ", fitkappa=", round(fitkappa, 3),
                                                   ", Loglik=", round(logLik(fit), 3)))

######################################################################
## Mpox
install.packages("readxl")
install.packages("lubridate")

library(readxl)
library(lubridate)

data <- read_excel("Desktop/Uni/BIO_4C12/SupplementaryFile1.xlsx", sheet = "3a. Incubation period_Pre-2022")
ExposureL <- ymd(data$ExposureL)
ExposureR <- ymd(data$ExposureR)

avgExp <- ifelse(!is.na(ExposureL) & !is.na(ExposureR), ExposureL + days(as.integer((ExposureR-ExposureL)/2)), 
                 ifelse(!is.na(ExposureL), ExposureL, ExposureR))
avgExp <- as.Date(avgExp, origin = "1970-01-01")
######################################################################
# rabies
load("/Users/ningruixie/Desktop/Uni/BIO_4C12/intervals.rda")
interval_merge <- as.data.frame(interval_merge)
SI <- subset(interval_merge, Type=="Serial Interval" & Days <100)

# why can't us use direct calculated mean & kappa?
interval = SI$Days
startPar = list(logmean = 3, logkappa = -0.1)
fit <- mle2(gPE.nll, 
            data = list(interval = interval),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = pe_interval,
                 perlang = dperlang(time, mean, kappa),
                 fit_perlang = dperlang(time, fitmean, fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=perlang, color="PErlang"), linewidth=1.5) +
  geom_line(aes(x=Time, y=fit_perlang, color="PErlang after fit"), linewidth=1.5) +
  labs(x = "Interval", y = "Count", title = paste0("PErlang -> PErlang |", 
                                                   " fitmean=", round(fitmean, 3), 
                                                   ", fitkappa=", round(fitkappa, 3),
                                                   ", Loglik=", round(logLik(fit), 3)))




############### check
kappa=0.3
mean=7
t <- 100
ts <- 0.05
time <- timeSeq(ts, t)
df <- data.frame(Time = time,
                 gamma = dgamma(time, shape=1/kappa, scale=mean*kappa),
                 perlang = dperlang(time, mean, kappa))
ggplot(df, aes(x=Time)) +
  geom_line(aes(y=gamma, color='Gamma'), linewidth=1.5, linetype='dashed') +
  geom_line(aes(y=perlang, color='Pseudo Erlang'), linewidth=1) +
  labs(title="Check: same mean & kappa")

ggplot(df, aes(x=Time)) +
  geom_line(aes(y=log(gamma), color='Gamma(l)'), linewidth=1.5, linetype='dashed') +
  geom_line(aes(y=log(perlang), color='Pseudo Erlang(l)'), linewidth=1) +
  labs(title="Check: same mean & kappa")


## area under the curve
sum(df$gamma*ts)
sum(df$perlang*ts)

quit()
num_mean <- sum(df$Time * df$perlang * 1) # ts = 1
num_var <- sum(df$Time^2 * df$perlang * 1) - num_mean^2
num_kappa <- num_var / num_mean^2
num_kappa

temp <- dperlang(time, 7, 0.02)
temp2 <- dgamma(time, shape=1/0.02, scale=7*0.02)

temp_df <- data.frame(Time=time, gamma=temp2, perlang=temp)
ggplot(temp_df, aes(x=Time)) +
  geom_line(aes(y=gamma,color='gamma')) +
  geom_line(aes(y=perlang,color='perlang'))

tempmean <- sum(temp2 * time * 0.012776)
tempvar <- sum(time^2 * temp2 * 0.012776) - tempmean^2
tempkappa <- tempvar/tempmean^2
tempkappa








