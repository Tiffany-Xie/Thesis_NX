library(deSolve)
library(devtools)
install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(ggplot2)#; theme_set(theme_bw(base_size=14))

######################################################################

nPE <- 12
ts <- 0.01
T <- 200
mu <- 10
μ <- 0

###################################################################### GOOD

tbound <- seq(0, T, by=ts)

plotSeries <- function(gamm, ode, tbound, title="Comparison plot", tmax=NA) {
  inc <- data.frame(
    Time = midpoints(tbound)
    , Gamma=diff(gamm)
    , ODE = diff(ode)
  )
  return(ggplot(inc)
         + aes(x=Time) 
         + geom_line(aes(y=Gamma, color = "Gamma"))
         + geom_line(aes(y=ODE, color = "ODE"))
         + xlim(c(0, tmax))
         + labs(y = "Incidence", title=title)
  )
} 


odeComp <- function(mu, kappa, tbound, nPE=12, plotMax=NA, pc=TRUE){
  
  ode <- Integration(nPE, mu, kappa, ts, T)
  gamm <- pgamma(tbound, 1/kappa, 1/(mu*kappa))
  
  if(pc){
    print(parCheck(gamm, ts, T))
    print(parCheck(ode, ts, T))
  }
  
  return(plotSeries(ode, gamm, tbound, tmax=plotMax))
}

twoPlots <- function(mu, kappa, tbound, plotMax=50, ymin=1e-6){
  oc <- (odeComp(mu, kappa, tbound=tbound, plotMax=plotMax)
         + ggtitle(paste("kappa = ", kappa))
  )
  print(oc)
  print(oc + scale_y_log10(limits=c(ymin, NA)))
}

twoPlots(10, 1/8, tbound)
#twoPlots(10, 1/4, tbound)
#twoPlots(10, 1/2, tbound)
#twoPlots(10, 3/4, tbound)
#twoPlots(10, 9/10, tbound)
#twoPlots(10, 0.99, tbound)

######################################################################

SInR <- function(time, states, params){
  with(as.list(c(params)), {
    S <- states[[1]]
    I <- states[2:(n+1)]
    R <- states[[n+2]]
    
    dS <- μ - β*S*sum(I) - μ*S
    dI <- c()
    dI[1] <- β*S*sum(I)-(n*γ + μ)*I[[1]]
    Iprev <- c(I[1:(n-1)])
    dI[2:n] <- n*γ*Iprev - (n*γ + μ)*I[2:n]
    dR <- n*γ*I[[n]] - μ*R
    return(list(c(dS, dI, dR)))
  })
}

n = 4
time <- seq(0,100,by=0.1)
params <- c(β=0.1, γ = 1/10, μ = 0.001, n = n)
states <- c(0.9, 0.1, numeric(n-1), 0)
names(states) <- c("S", paste0("I", 1:n), "R") 

soln <- ode(y = states,
            times = time,
            func = SInR,
            parms = params)
#plot(soln)
plot(time, soln[,"R"], type = "l")
head(rowSums(soln[,-1]))  # Checking

######################################################################

sinner <- function(time, states, params) {
  with(as.list(c(params)), {
    S <- states[[1]]
    I <- states[2:(n+1)]
    R <- states[[n+2]]
    
    dS <- μ - β*S*sum(I) - μ*S
    dI <- c()
    dI[1] <- β*S*sum(I)-(a + μ)*I[[1]]
    Iprev <- c(I[1:(n-1)])
    dI[2:n] <- a*r^(0:(n-2))*Iprev - (a*r^(1:(n-1)) + μ)*I[2:n]
    dR <- a*r^(n-1)*I[[n]] - μ*R
    return(list(c(dS, dI, dR)))
  })
}

#in flow and out flow of the model separate.
######################################################################
######################################################################
angel <- function(){
  rates <- something
  flows <-  something
  flowPrev <-  something
  dI <- something_about_the_flows
  dI[[1]] <- dI[[1]] + extraTerm
}

n = 12
mu = 10
kappa = 1/4
r <- kappa2r(kappa, n)
a <- (1-1/r^n)/(mu*(1-1/r))

states <- c(0.9, 0.1, numeric(n-1), 0)
names(states) <- c("S", paste0("I", 1:n), "R") 
params <- c(β=0.1, μ = 0.001, a = a, r = r, n=n)
soln <- ode(y = states,
            times = time,
            func = sinner,
            parms = params)
#plot(soln)
plot(time, soln[,"R"], type = "l")
head(rowSums(soln[,-1]))

######################################################################

SEmInR <- function(time, states, params){
  with(as.list(c(params)), {
    S <- states[[1]]
    E <- states[2:(m+1)]
    I <- states[(m+2):(m+n+1)]
    R <- states[[m+n+2]]
    
    dS <- μ - β*S*sum(I) - μ*S
    dE <- c()
    dE[1] <- β*S*sum(I)-(m*σ + μ)*E[[1]]
    Eprev <- c(E[1:(m-1)])
    dE[2:m] <- m*σ*Eprev - (m*σ + μ)*E[2:m]
    
    dI <- c()
    dI[1] <- m*σ*E[m] - (n*γ + μ)*I[[1]]
    Iprev <- c(I[1:(n-1)])
    dI[2:n] <- n*γ*Iprev - (n*γ + μ)*I[2:n]
    
    dR <- n*γ*I[[n]] - μ*R
    return(list(c(dS, dE, dI, dR)))
  })
}

m = 4
n = 3
params <- c(β=0.1, γ = 1/10, σ = 0.2 ,μ = 0.001, n = n, m = m)
states <- c(0.9, numeric(m), 0.1, numeric(n-1), 0)
names(states) <- c("S", paste0("E", 1:m), paste0("I", 1:n), "R") 

soln <- ode(y = states,
            times = time,
            func = SEmInR,
            parms = params)
#plot(soln)
plot(time, soln[,"S"], type = "l")
head(rowSums(soln[,-1])) # Checking

######################################################################
seminar <- function(time, states, params) {
  with(as.list(c(params)), {
    S <- states[[1]]
    E <- states[2:(m+1)]
    I <- states[(m+2):(m+n+1)]
    R <- states[[m+n+2]]
    
    dS <- μ - β*S*sum(I) - μ*S
    dE <- c()
    dE[1] <- β*S*sum(I)-(a1 + μ)*E[[1]]
    Eprev <- c(E[1:(m-1)])
    dE[2:m] <- a1*r1^(0:(m-2))*Eprev - (a1*r1^(1:(m-1)) + μ)*E[2:m]
    
    dI <- c()
    dI[1] <- a1*r1^(m-1)*E[m] - (a2 + μ)*I[[1]]
    Iprev <- c(I[1:(n-1)])
    dI[2:n] <- a2*r2^(0:(n-2))*Iprev - (a2*r2^(1:(n-1)) + μ)*I[2:n]
    
    dR <- a2*r2^(n-1)*I[[n]] - μ*R
    return(list(c(dS, dE, dI, dR)))
  })
}

m=6
n=7
mu1=5 # sigma
mu2 = 10 # gamma
kappa1 = 1/4
kappa2 = 1/3
r1 <- kappa2r(kappa1, m)
a1 <- (1-1/r1^m)/(mu1*(1-1/r1))
r2 <- kappa2r(kappa2, n)
a2 <- (1-1/r2^n)/(mu2*(1-1/r2))

params <- c(β=0.1,μ = 0.001, a1=a1, a2=a2, r1=r1, r2=r2, n = n, m = m)
states <- c(0.9, numeric(m), 0.1, numeric(n-1), 0)
names(states) <- c("S", paste0("E", 1:m), paste0("I", 1:n), "R") 

soln <- ode(y = states,
            times = time,
            func = seminar,
            parms = params)
#plot(soln)
plot(time, soln[,"S"], type = "l")
head(rowSums(soln[,-1])) # Checking


######################################################################
cgs <- function(r, n){
  delta <- (1:n) - (n+1)/2
  return(exp(delta*log(r)))
}

geomFlow <- function(mu, n, kappa, ts, T){
  r <- kappa2r(kappa, n)
  D = cgs(r, n)
  D = D*mu/(sum(D))
  return(1/D)
}

ts <- 0.1
n <- 4
T <- 200
mu <- 5
pen <- 12
kappa <- 1/4

geomFlow(mu, pen, kappa, ts, T)



 