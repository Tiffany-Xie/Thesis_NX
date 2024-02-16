library(deSolve)
library(pseudoErlang)
library(gridExtra)
library(ggplot2) 
theme_set(theme_minimal(base_size = 17))

###################################################################### Figure1 - ErlangVsOde

erlang <- function(x, n, γ) {
  (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

InR <- function(t, states, params) {
  with(as.list(c(params)), {
    I <- states[1:n]
    R <- states[[n+1]]
    
    Iprev <- c(0, I[1:(n-1)])
    dI <- n*γ*Iprev - (n*γ + μ)*I
    dR <- n*γ*I[[n]] - μ*R
    return(list(c(dI, dR)))
  })
}

compare <- function(n, γ, μ=0, ts, T) {
  time <- timeSeq(ts, T, FALSE)
  timem <- timeSeq(ts, T)
  
  params <- c(n=n, γ=γ, μ=0)
  states <- c(1, numeric(n))
  names(states) <- c(paste0("I", 1:n), "R")
  
  erlang <- erlang(timem, n, γ)
  ode <- ode(y = states,
             times = time,
             func = InR,
             parms = params)
  
  df <- data.frame(Time = timem,
                   Erlang = erlang,
                   ODE = diff(ode[,"R"])/ts)
  
  ggplot(df, aes(x=Time)) +
    geom_line(aes(y = Erlang, col = "Erlang"), linewidth = 1) +
    geom_line(aes(y = ODE, col = "ODE"), linewidth = 2, linetype = "twodash") +
    geom_vline(xintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
    labs(x = "Stage duration (days)",
         y = "Probability Density",
         title = paste("n = ", n))
}


ts <- 0.01
T <- 40
γ <- 0.1
N <- c(4,8,12,18) # 4 8 12 18

plots <- lapply(N, function(N) compare(N, γ, μ=0, ts, T))
grid.arrange(grobs = plots, ncol = 2, nrow = 2)


###################################################################### Figure2 - PossibleMatch

InR_geom <- function(t, states, params) {
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
  df <- expand.grid(Time = (time[-1]+time[-length(time)])/2, nPE = nPE, γ = γ, r = r)
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
    PPE <- c(PPE, diff(soln[,"R"])/diff(time))
  }
  df$PPE <- PPE
  return(df)
}


Cfdens <- function(dfE, dfPE) {
  ggplot(dfPE, aes(x=Time, y = PPE, color=factor(r))) + 
    geom_line() +
    scale_color_manual(values = viridis(length(unique(r)))) +
    geom_vline(xintercept = 1/dfPE$γ[1], color = "red", linetype = "dashed", linewidth = 1) +
    geom_line(data=dfE, aes(x = Time, y = PE), color = "black", linewidth = 1.5) +
    labs(x = "Stage duration (days)",
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
    geom_hline(yintercept = 1/(nE), color = "red", linetype = "dashed", linewidth = 1)
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
    labs(y="mean")
  
  p2 <- ggplot(dfPEK, aes(x=r, y=PEK)) + geom_line() +
    geom_hline(yintercept=EK, color="red", linetype = "dashed", linewidth = 1) +
    labs(y="K")
  grid.arrange(p1, p2, ncol=2)
}


# Simulation ####
library(viridis)
time = seq(0,40,by=0.01) 
γ <- 0.1 
μ <- 0
nPE <- 12
ts <- 0.01

# nE == 2
r <- seq(2.5, 3.5, by=0.1)
nE <- 2
dfE <- Edens(time, γ, nE)
dfPE <- PEdens(time, γ, μ, r, nPE, InR_geom)
dfK <- Kappaf(r, γ, nPE)
p1 <- Cfdens(dfE, dfPE) + ggtitle(paste("n = ", nE))
p2 <- Cfkappa(dfK) + ggtitle("")
grid.arrange(p1, p2, ncol = 2) # r (2.75, 3.25)
# Num
#Cfnum(dfE, dfPE, r, ts) # not the correct one for numerical calculation

# nE == 4
r <- seq(1.2, 2.2, by=0.1)
nE <- 4
dfE <- Edens(time, γ, nE)
dfPE <- PEdens(time, γ, μ, r, nPE, InR_geom)
dfK <- Kappaf(r, γ, nPE)
p1 <- Cfdens(dfE, dfPE) + ggtitle(paste("n = ", nE))
p2 <- Cfkappa(dfK) + ggtitle("")
grid.arrange(p1, p2, ncol = 2) #r (1.5, 1.75)
# Num
#Cfnum(dfE, dfPE, r, ts)

# nE == 8
r <- seq(1.2, 1.30, by=0.01)
nE <- 8
dfE <- Edens(time, γ, nE)
dfPE <- PEdens(time, γ, μ, r, nPE, InR_geom)
dfK <- Kappaf(r, γ, nPE)
p1 <- Cfdens(dfE, dfPE) + ggtitle(paste("n = ", nE))
p2 <- Cfkappa(dfK) + ggtitle("")
grid.arrange(p1, p2, ncol = 2) # r (1.2, 1.25)
# Num
#Cfnum(dfE, dfPE, r, ts)

###################################################################### Figure3 TheoreticalAlignment

compare <- function(nfix, mu, n, ts, T) {
  time <- timeSeq(ts,T,FALSE)
  timem <- timeSeq(ts,T)
  ode <- Integration(nfix, mu, 1/n, ts, T)
  gamm <- pgamma(time, n, n/mu)
  
  df <- data.frame(Time = timem,
                   ODE = diff(ode)/ts,
                   Gamma = diff(gamm)/ts)
  ggplot(df, aes(x=Time)) + 
    geom_line(aes(y=ODE, color = "ODE"), linewidth=1) +
    geom_line(aes(y=Gamma, color = "Erlang"), linewidth=1) + 
    geom_vline(xintercept = mu, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = paste("n = ", n),
         x = "Stage Duration, (days)",
         y = "Probability Density")
}

nfix <- 12
mu <- 10
n <- c(2, 4, 8, 10)
ts <- 0.01
T <- 70

plots <- lapply(n, function(n) compare(nfix, mu, n, ts, T))
grid.arrange(grobs = plots, ncol = 2, nrow = 2)

###################################################################### Figure4 SinrVsSigr
library(tidyr)
library(dplyr)

Integration <- function(params, states, ts, T, model) {
  time <- timeSeq(ts, T, FALSE)
  soln <- ode(y = states,
              times = time,
              func = model,
              parms = params)
  return(soln)
}

compare <- function(flat, geom, ts, T, title) {
  time <- timeSeq(ts, T, FALSE)
  df <- data.frame(Time = time, 
                   FlatI = rowSums(flat[,3:(ncol(flat)-1)]),
                   GeomI = rowSums(geom[,3:(ncol(geom)-1)]))
  ggplot(df, aes(x=Time)) +
    geom_line(aes(y=FlatI, color = "Flat")) +
    geom_line(aes(y=GeomI, color = "Geom")) +
    labs(y = "Density", 
         title = paste(title, "Infectious"))
}

##############################################

SIR <- function(time, states, params) {
  with(params, {
    stopifnot(length(outrate)==n)
    S <- states[[1]]
    I <- states[2:(n+1)]
    R <- states[[n+2]]
    
    outflow <- outrate*I
    inrate <- outrate-μ
    inflow <- c(β*S*sum(I), (I*inrate)[1:(n-1)])
    
    dS <- μ - β*S*sum(I) - μ*S
    dI <- inflow - outflow
    dR <- inrate[n]*I[[n]] - μ*R
    return(list(c(dS, dI, dR)))
  })
}

SEIR <- function(time, states, params) {
  with(params, {
    stopifnot(length(outrate)==n+m)
    S <- states[[1]]
    #E <- states[2:(m+1)]
    #I <- states[(m+2):(m+n+1)]
    mid <- states[2:(m+n+1)]
    R <- states[[m+n+2]]
    
    outflow <- outrate*mid
    inrate <- outrate-μ
    inflow <- c(β*S*sum(mid[(m+1):(m+n)]), (mid*inrate)[1:(m+n-1)])
    
    dS <- μ - β*S*sum(mid[(m+1):(m+n)]) - μ*S
    dE <- (inflow - outflow)[1:m]
    dI <- (inflow - outflow)[(m+1):(m+n)]
    dR <- inrate[m+n]*mid[[m+n]] - μ*R
    return(list(c(dS, dE, dI, dR)))
  })
}

#############################################

SInRFlow <- function(β, mu, n, μ, ts, T) {
  gamma <- 1/mu
  outrate <- rep(n*gamma + μ, times=n)
  
  params <- list(β=β, n=n, μ=μ, outrate=outrate)
  states <- c(0.9, 0.1, numeric(n-1), 0)
  names(states) <- c("S", paste0("I", 1:n), "R")
  return(Integration(params, states, ts, T, SIR))
}

sinnerFlow <- function(β, mu, kappa, n, μ, ts, T) {
  r <- kappa2r(kappa, n)
  a <- (1-1/r^n)/(mu*(1-1/r))
  print(c(r, a))
  outrate <- a*r^(0:(n-1)) + μ
  
  params <- list(β=β, n=n, μ=μ, outrate=outrate)
  states <- c(0.9, 0.1, numeric(n-1), 0)
  names(states) <- c("S", paste0("I", 1:n), "R")
  return(Integration(params, states, ts, T, SIR))
}

SEmInRFlow <- function(β, muE, muI, m, n, μ, ts, T) {
  gammaE <- 1/muE
  gammaI <- 1/muI
  outrate <- c(rep(m*gammaE + μ, times=m),
               rep(n*gammaI + μ, times=n))
  
  params <- list(β=β, m=m, n=n, μ=μ, outrate=outrate)
  states <- c(0.9, numeric(m), 0.1, numeric(n-1), 0)
  names(states) <- c("S", paste0("E", 1:m), paste0("I", 1:n), "R")
  return(Integration(params, states, ts, T, SEIR))
}

seminarFlow <- function(β, muE, muI, kappaE, kappaI, m, n, μ, ts, T) {
  rE <- kappa2r(kappaE, m)
  aE <- (1-1/rE^m)/(muE*(1-1/rE))
  rI <- kappa2r(kappaI, n)
  aI <- (1-1/rI^n)/(muI*(1-1/rI))
  outrate <- c(aE*rE^(0:(m-1)) + μ, aI*rI^(0:(n-1)) + μ)
  
  params <- list(β=0.1,μ = 0.001, n = n, m = m, outrate = outrate)
  states <- c(0.9, numeric(m), 0.1, numeric(n-1), 0)
  names(states) <- c("S", paste0("E", 1:m), paste0("I", 1:n), "R")  
  return(Integration(params, states, ts, T, SEIR))
}

###############################################

β <- 0.1
mu <- 10
kappa <- 1/4
fixn <- 12
μ <- 0.001
ts <- 0.1
T <- 70

n <- 2
kappa <- 1/2
time <- timeSeq(ts, T, FALSE)
sinner <- sinnerFlow(β, mu, kappa, fixn, μ, ts, T)
sinr <- SInRFlow(β, mu, n, μ, ts, T)
dfn2 <- data.frame(Time = time, 
                   Flat = rowSums(sinr[,3:(ncol(sinr)-1)]),
                   Geom = rowSums(sinner[,3:(ncol(sinner)-1)]))
n <- 4
kappa <- 1/4
time <- timeSeq(ts, T, FALSE)
sinner <- sinnerFlow(β, mu, kappa, fixn, μ, ts, T)
sinr <- SInRFlow(β, mu, n, μ, ts, T)
dfn4 <- data.frame(Time = time, 
                   Flat = rowSums(sinr[,3:(ncol(sinr)-1)]),
                   Geom = rowSums(sinner[,3:(ncol(sinner)-1)]))

dfn2$id <- 'n=2'
dfn4$id <- 'n=4'
combined_df <- rbind(dfn2, dfn4)
long_df <- combined_df %>% 
  pivot_longer(cols = c(Flat, Geom), names_to = "variable", values_to = "value")


ggplot(long_df, aes(x = Time, y = value, color = id, linetype = variable)) +
  geom_line(linewidth=1) +
  scale_color_manual(values = c('n=2' = 'blue', 'n=4' = 'orange')) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(color = "n", linetype = "Model Type",
       x = "Time, (days)",
       y = "Density")


###################################################################### Figure5 - VersatileK

nfix <- 12
mu <- 10
ts <- 0.01
T <- 60
kappa <- c(0.25, 0.45, 0.65, 0.85)
time <- timeSeq(ts, T)

df <- expand.grid(Time = time, Kappa = kappa)

Ode <- c()
for (k in kappa) {
  ode <- diff(Integration(nfix, mu, k, ts, T))/ts
  Ode <- c(Ode, ode)
}

df[["ODE"]] <- Ode

ggplot(df, aes(x=Time, y=ODE, color = factor(Kappa))) + geom_line(linewidth = 1) +
  geom_vline(xintercept = mu, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Stage duration (days)",
       y = "Probability Density")






###################################################################### 


  