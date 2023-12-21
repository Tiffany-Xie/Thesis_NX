library(deSolve)
library(devtools)
#install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(gridExtra)
library(ggplot2)#; theme_set(theme_bw(base_size=14))

###################################################################### 

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
                   FlatS = flat[,"S"],
                   GeomS = geom[,"S"],
                   FlatR = flat[,"R"],
                   GeomR = geom[,"R"])
  p1<- ggplot(df, aes(x=Time)) +
    geom_line(aes(y=FlatS, color = "Flat")) +
    geom_line(aes(y=GeomS, color = "Geom")) +
    labs(y = "Density", 
         title = paste(title, "(S)"))
  p2<- ggplot(df, aes(x=Time)) +
    geom_line(aes(y=FlatR, color = "Flat")) +
    geom_line(aes(y=GeomR, color = "Geom")) +
    labs(y = "Density", 
         title = paste(title, "(R)"))
  grid.arrange(p1, p2, ncol = 2)
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

######################################################################

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

######################################################################

β <- 0.1
mu <- 10
kappa <- 1/4
fixn <- 12
n <- 4
μ <- 0.001
ts <- 0.1
T <- 100

sinner <- sinnerFlow(β, mu, kappa, fixn, μ, ts, T)
#head(rowSums(soln[,-1]))
sinr <- SInRFlow(β, mu, n, μ, ts, T)

compare(sinr, sinner, ts, T, title="SInR & Geom SInR")

######################################################################

β <- 0.1
muE <- 5
muI <- 10
kappaE <- 1/4
kappaI <- 1/3
fixm <- 6
fixn <- 7
m <- 4
n <- 3
μ <- 0.001

seminar <- seminarFlow(β, muE, muI, kappaE, kappaI, fixm, fixn, μ, ts, T)
#head(rowSums(soln[,-1]))
seminr <- SEmInRFlow(β, muE, muI, m, n, μ, ts, T)

compare(seminr, seminar, ts, T, title="SEmInR & Geom SEmInR")




