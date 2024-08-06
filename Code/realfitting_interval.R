library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)
library(pseudoErlang)

library(shellpipes)
loadEnvironments()

######################################################################
## rabies
######################################################################
url <- "https://github.com/eliminaterabies/R0paper/raw/main/public_data/intervals.rda"

loadData <- function(url) {
  temp <- tempfile()
  download.file(url, temp, mode = "wb")
  load(temp, envir = .GlobalEnv)
  unlink(temp)
}
loadData(url)

# visualization
######################################################################

interval_merge <- as.data.frame(interval_merge)
SI <- subset(interval_merge, Type=="Serial Interval")
time <- seq(0, max(SI$Days), length.out=dim(SI)[1])
SI$Time <- time
SI$theogamma <- dgamma(time, shape=1/1.935, scale=1.935* 28.2)

ggplot(SI) +
  geom_histogram(aes(x=Days, y=after_stat(density))) +
  labs(title='Serial Interval (density) | totMean = 28.2; totKappa = 1.935')

SI <- subset(interval_merge, Type=="Serial Interval" & Days <100)
time <- seq(0, max(SI$Days), length.out=dim(SI)[1])
SI$Time <- time
SI$theogamma <- dgamma(time, shape=1/0.696, scale=0.696*21.56)
SI$theoperlang <- dperlang(time, 21.56, 0.696)

ggplot(SI) +
  geom_histogram(aes(x=Days, y=after_stat(density))) +
  geom_line(aes(x=Time, y=theogamma), linewidth=1) +
  labs(title='Serial Interval < 100 & Gamma| totMean = 21.56; totKappa = 0.696')

ggplot(SI) +
  geom_histogram(aes(x=Days, y=after_stat(density))) +
  geom_line(aes(x=Time, y=theoperlang), linewidth=1) +
  labs(title='Serial Interval < 100 & PErlang | totMean = 21.56; totKappa = 0.696')

# Serial Interval: gamma fit 
######################################################################

SI$Days[SI$Days==0] <- 0.01

startPar = list(logmean = 4, logkappa = -0.4)
fit <- mle2(gg.nll, 
            data = list(interval = SI$Days),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = SI$Days,
                 fit_gamma = dgamma(time, shape=1/fitkappa, scale=fitmean*fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=fit_gamma), linewidth=1) +
  labs(x = "Interval", y = "Density", title = paste0("Serial Interval -> Gamma |", 
                                                   " fitmean=", round(fitmean, 3), 
                                                   ", fitkappa=", round(fitkappa, 3)))

## why can't us use direct calculated mean & kappa?

# Serial Interval: pseudo erlang fit 
######################################################################

startPar = list(logmean = 3, logkappa = -1.1)
fit <- mle2(gPE.nll, 
            data = list(interval = SI$Days),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(exp(coef(fit)[["logmean"]]))
fitkappa = print(exp(coef(fit)[["logkappa"]]))
df <- data.frame(Time = time,
                 interval = SI$Days,
                 fit_perlang = dperlang(time, fitmean, fitkappa))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=fit_perlang), linewidth=1) +
  labs(x = "Interal, days", y = "Density", title = paste0("Serial Interval -> PErlang |", 
                                                   " fitmean=", round(fitmean, 3), 
                                                   ", fitkappa=", round(fitkappa, 3)))

# Serial Interval: lognorm fit
######################################################################
# SI <- subset(interval_merge, Type=="Serial Interval")
# time <- seq(0, max(SI$Days), length.out=dim(SI)[1])
# SI$Time <- time
# SI$Days[SI$Days==0] <- 0.01

gln.nll <- function(mean, sd, interval) {
  -sum(dlnorm(interval, meanlog=log(mean^2/sqrt(mean^2 + sd^2)), sdlog=log(1+sd^2/mean^2), log=TRUE))
}

startPar = list(mean = 21, sd = 20)
fit <- mle2(gln.nll, 
            data = list(interval = SI$Days),
            start = startPar,
            method = "Nelder-Mead",
            control = list(maxit = 10000))

fitmean = print(coef(fit)[["mean"]])
fitsd = print(coef(fit)[["sd"]])
df <- data.frame(Time = time,
                 interval = SI$Days,
                 fit_lnorm = dlnorm(time, meanlog=log(fitmean^2/sqrt(fitmean^2 + fitsd^2)), sdlog=log(1+fitsd^2/fitmean^2)))
ggplot(df) + 
  geom_histogram(aes(x=interval, y = after_stat(density))) +
  geom_line(aes(x=Time, y=fit_lnorm), linewidth=1) +
  labs(x = "Interal, days", y = "Density", title = paste0("Serial Interval -> LogNorm |", 
                                                        " fitmean=", round(fitmean, 3), 
                                                        ", fitsd=", round(fitsd, 3)))
## how can we evaluate the fitting result 

quit()
######################################################################
## Mpox: + Triangular distribution
######################################################################

library(readxl)
library(lubridate)

data <- read_excel("Desktop/Uni/BIO_4C12/SupplementaryFile1.xlsx", sheet = "3a. Incubation period_Pre-2022")
ExposureL <- ymd(data$ExposureL)
ExposureR <- ymd(data$ExposureR)

avgExp <- ifelse(!is.na(ExposureL) & !is.na(ExposureR), ExposureL + days(as.integer((ExposureR-ExposureL)/2)), 
                 ifelse(!is.na(ExposureL), ExposureL, ExposureR))
avgExp <- as.Date(avgExp, origin = "1970-01-01")


# -sum(Integral(dgamma * trangularPDF))
# dgamma(x); trangular(x, a, b, c)

gPE.intv.nll <- function(logmean, logkappa, interval) { # interval=c(l,r)
  left <- interval[1]
  right <- interval[2]
  f <- function(x) TxPE(x, exp(mean), exp(kappa), left, right)
  -sum(log(integrate(f, lower=left, upper=right)))
}

dtriangular <- function(x, left, right, mode=0) {
  if (mode==0) mode <- (left+right)/2
  if (x < left || x > right) return(0)
  if (x <= mode) return((2 * (x - left)) / ((right - left) * (mode - left)))
  return((2 * (right - x)) / ((right - left) * (right - mode)))
}

TxPE <- function(x, mean, kappa, left, right,) {
  dtriangular(x, left, right) * dperlang(x, mean, kappa)
}






