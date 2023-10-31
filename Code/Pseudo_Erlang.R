# Packages ####
library(viridis)
library(ggplot2)
library(deSolve)
library(gridExtra)

library(shellpipes)
startGraphics(width=11)

# Models ####
erlang <- function(x, n, γ) {
  (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

SInR_geom <- function(t, states, params) {
  with(as.list(c(params)), {
    I <- states[1:n]
    R <- states[[n+1]]
    
    Iprev <- c(0, I[1:(n-1)])
    dI <- a*r^(c(0,0:(n-2)))*Iprev - (a*r^(0:(n-1))+μ)*I
    dR <- a*r^(n-1)*I[[n]] - μ*R
    return(list(c(dI, dR)))
  })
}

# Functions ####
Edens <- function(time, γ, nE) {
  df <- data.frame(Time = time)
  df$PE <- erlang(time, nE, γ)
  return(df)
}

PEdens <- function(time, γ, μ, r, nPE, model) {
  df <- expand.grid(Time = time, nPE = nPE, γ = γ, r = r)
  df$a <- with(df, (1/r^nPE-1)*γ/(1/r-1))
  states <- c(1, numeric(nPE))
  names(states) <- c(paste0("I", 1:nPE),"R")
  PPE <- c()
  for (i in r) {
    a <- (1/i^nPE-1)*γ/(1/i-1)
    params <- c(γ = γ, μ = μ, n = nPE, a = a, r = i)
    soln <- ode(y = states,
                times = time, 
                func = model, 
                parms = params)
    PPE <- c(PPE, c(diff(soln[,"R"])/diff(time), NA))
  }
  df$PPE <- PPE
  return(na.omit(df)) ## Is this a trap for future Ningrui?
}


Cfdens <- function(dfE, dfPE) {
  ggplot(dfPE, aes(x=Time, y = PPE, color=factor(r))) + 
    geom_line() +
    scale_color_manual(values = viridis(length(unique(r)))) +
    geom_vline(xintercept = 1/dfPE$γ[1], color = "red", linetype = "dashed", linewidth = 1) +
    geom_line(data=dfE, aes(x = Time, y = PE), color = "black", linewidth = 1.5) +
    labs(title = paste0("Compare different r: D=", 1/dfPE$γ[1], ", nE=", nE, ", nPE=", nPE),
         x = "Stage duration (days)",
         y = "Probability Density")
}


Kappaf <- function(r, γ, nPE) {
  df <- data.frame(r = r)
  df$a <- with(df, (1/r^nPE-1)*γ/(1/r-1))
  df$mean <- with(df, 1/a*(1/r^nPE-1)/(1/r-1))
  df$var <- with(df, 1/a^2*(1/r^(2*nPE)-1)/(1/r^2-1))
  df$K <- with(df, var/mean^2)
  return(df)
}

Cfkappa <- function(dfK) {
  ggplot(dfK, aes(x=r, y=K)) + geom_line() +
    geom_hline(yintercept = 1/(nE), color = "red", linetype = "dashed", linewidth = 1) +
    labs(title="Coefficient of Variation^2 (K)")
}


Cfnum <- function(dfE, dfPE, r, ts) {
  Emean = sum(dfE$Time * dfE$PE * ts)
  Evar = sum(dfE$Time^2 * dfE$PE * ts) - Emean^2
  EK = Evar/Emean^2
  dfPEmean <- data.frame(r = r)
  dfPEK <- data.frame(r = r)
  tempm <- c()
  tempK <- c()
  for (i in r) {
    rows <- dfPE[dfPE$r == i,]
    PEmean <- sum(rows$Time * rows$PPE * ts)
    PEvar <- sum(rows$Time^2 * rows$PPE * ts) - PEmean^2
    PEK = PEvar/PEmean^2
    tempm <- c(tempm, PEmean)
    tempK <- c(tempK, PEK)
  }
  dfPEmean$PEmean <- tempm
  dfPEK$PEK <- tempK
  
  p1 <- ggplot(dfPEmean, aes(x=r, y=PEmean)) + geom_line() +
    geom_hline(yintercept=Emean, color="red", linetype = "dashed", linewidth = 1) +
    labs(title=paste0("Compare numerically calculated mean, nE=", nE, ", nPE=", nPE),
         y="mean")
  
  p2 <- ggplot(dfPEK, aes(x=r, y=PEK)) + geom_line() +
    geom_hline(yintercept=EK, color="red", linetype = "dashed", linewidth = 1) +
    labs(title=paste0("Compare numerically calculated K, nE=", nE, ", nPE=", nPE),
         y="K")
  grid.arrange(p1, p2, ncol=2)
}


# Simulation ####
time = seq(0,40,by=0.1) 
γ <- 0.1 
μ <- 0
nPE <- 12
ts <- 0.1

# nE == 2
r <- seq(2.5, 3.5, by=ts)
nE <- 2
dfE <- Edens(time, γ, nE)
dfPE <- PEdens(time, γ, μ, r, nPE, SInR_geom)
dfK <- Kappaf(r, γ, nPE)
p1 <- Cfdens(dfE, dfPE)
p2 <- Cfkappa(dfK)
grid.arrange(p1, p2, ncol = 2) # r (2.75, 3.25)
# Num
Cfnum(dfE, dfPE, r, ts)

# nE == 4
r <- seq(1.5, 2.5, by=0.1)
nE <- 4
dfE <- Edens(time, γ, nE)
dfPE <- PEdens(time, γ, μ, r, nPE, SInR_geom)
dfK <- Kappaf(r, γ, nPE)
p1 <- Cfdens(dfE, dfPE)
p2 <- Cfkappa(dfK)
grid.arrange(p1, p2, ncol = 2) #r (1.5, 1.75)
# Num
Cfnum(dfE, dfPE, r, ts)

# nE == 8
r <- seq(1.2, 1.30, by=0.01)
nE <- 8
dfE <- Edens(time, γ, nE)
dfPE <- PEdens(time, γ, μ, r, nPE, SInR_geom)
dfK <- Kappaf(r, γ, nPE)
p1 <- Cfdens(dfE, dfPE)
p2 <- Cfkappa(dfK)
grid.arrange(p1, p2, ncol = 2) # r (1.2, 1.25)
# Num
Cfnum(dfE, dfPE, r, ts)


# Calculate numerically from the output
#mean = sum(df$Time * df$P * ts)
#var = sum(df$Time^2 * df$P * ts) - mean^2
#K = var/mean^2

# Numerically find r to match κ=1/n_E
nPE <- 12
γ <- 0.1
nE <- 8

f <- function(r) 1/(γ*(1-1/r^nPE)/(1-1/r))^2 * (1-1/r^(2*nPE))/(1-1/r^2) * γ^2 - 1/nE
root <- uniroot(f, interval = c(1.2, 1.25))
r <- root$root # 1.241118

dfE <- Edens(time, γ, nE)
dfPE <- PEdens(time, γ, μ, r, nPE, SInR_geom)
dfK <- Kappaf(r, γ, nPE) # dfK$K ≈ 1/nE = 0.125
Cfdens(dfE, dfPE)

check <- sum(dfPE$PPE * ts)
dfPE$PPE = dfPE$PPE/check
mean <- sum(dfPE$Time * dfPE$PPE * ts) # 9.949945
var <- sum(dfPE$Time^2 * dfPE$PPE * ts) - mean^2 # 12.49957 
K = var/mean^2 # 0.1262565

print(c(check, mean, var, K))








