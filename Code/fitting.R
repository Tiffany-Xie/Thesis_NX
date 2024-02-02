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
  inc <- diff(sinr[,"inc"]) #/ts
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  df <- data.frame(inc = inc,
                   obs = obs)
  return(df)
}

fitting <- function(df, βe, De, n, μ, S0, I0, ts, T, seed, plot=TRUE) {
  obs = df$obs
  sir.nll <- function(βe, De, obs) {
    process_betae <<- c(process_betae, βe)
    process_De <<- c(process_De, De)
    #print(c(βe=βe, De=De))
    out <- as.data.frame(SInRFlow(β=exp(βe), D=exp(De), n=n, μ=μ, S0=S0, I0=I0, ts=ts, T=T))
    nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=1000, log=TRUE))
  }
  
  params0 <-list(βe=βe, De=De)
  set.seed(seed) 
  tryCatch({
    fit0 <- mle2(sir.nll, start=params0, data=list(obs=obs))
    if (plot) plotting(df, fit0, ts, T)
    return(fit0)
  }, error = function(e) {
    "That's unfortunate"
  }, finally = {
    1
  })
  
}

plotting <- function(df, fit, ts, T) {
  plot(timeSeq(ts, T), df$obs, type="l", xlab='Time, (Days)', ylab='I(t)')
  t <- timeSeq(ts, T)
  mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(coef(fit)[["βe"]]), D=exp(coef(fit)[["De"]]), n, μ, S0, I0, ts, T)))
  lines(diff(mod.prep$inc)~t, col = "red", lwd=3)
  lines(timeSeq(ts, T), df$inc, col = "blue")
}

## Always works ####################################################################     

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

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 77) # seed = 77
process_betae <- c()
process_De <- c()
ans <- fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, 76) # seed = 76


pardf <- data.frame(order = seq(0:(length(process_betae)-1)), 
                    betae = process_betae,
                    De = process_De)

p <- ggplot(pardf, aes(x = De, y = betae)) +
  geom_point(aes(color = order)) +  # Color by Time
  geom_path(alpha = 0.5) +  # Trace path
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient
  theme_minimal() +
  labs(x = "D", y = "Beta", color = "Time")

library(dplyr)
pardf$lead_D <- dplyr::lead(pardf$De)
pardf$lead_beta <- dplyr::lead(pardf$betae)

# Filter for every nth row to avoid cluttering, adjust n as needed
n <- 10  # For example, add an arrow every 10 points
arrow_data <- pardf[seq(1, nrow(pardf), by = n), ]

p + geom_segment(data = arrow_data, aes(xend = lead_D, yend = lead_beta), arrow = arrow(type = "closed", length = unit(0.15, "inches")), lineend = "round")

## Change time #################################################################### Change time

T <- 200

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 71) # seed = 71
fit = fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, 72) # seed = 72

## Change population (N) #################################################################### Change population (N)

T <- 100
S0 <- 99999

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 70) # seed = 70
fit = fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, 69) # seed = 69
 
## Change ts ####################################################################

S0 <- 999
ts <- 0.1 # works

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 67) # seed = 70
fit = fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, 66) # seed = 69

## Seeds mle #################################################################### Seeds - mle

seeds <- seq(1,100) 
ts <- 1

succ_mle <- c()
sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 67)

params0 <-list(βe=-0.5, De=2)

for (s in seeds) {
  ans <- tryCatch({
    fit <- fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, s, FALSE)
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

print(succ_mle)

## Seeds NB ####################################################################

seeds <- seq(1,100) 

succ_nb <- c()
sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)

params0 <-list(βe=-0.5, De=2)

for (s in seeds) {
  df = simObs(sinr, arp, nbs, s)
  
  ans <- tryCatch({
    fit <- fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, 33, FALSE) #mle seed = 33
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

print(succ_nb)


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



quit()
######################################################################







