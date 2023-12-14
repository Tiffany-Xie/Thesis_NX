library(deSolve)
library(devtools)
#install_github("Tiffany-Xie/pseudoErlang")
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
         + geom_line(aes(y=Gamma, color = "gamma"))
         + geom_line(aes(y=ODE, color = "ODE"))
         + xlim(c(0, tmax))
         + labs(y = "incidence", title=title)
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
twoPlots(10, 1/4, tbound)
twoPlots(10, 1/2, tbound)
twoPlots(10, 3/4, tbound)
twoPlots(10, 9/10, tbound)
twoPlots(10, 0.99, tbound)

######################################################################

SInR <- function(time, states, params){
  with(as.list(c(params)), {
    S <- states[[1]]
    I <- states[2:(n+1)]
    R <- states[[n+2]]
    
    dS <- μ - β*S*sum(I) - μ*S
    Iprev <- c(β*S*sum(I), I[1:(n-1)])
    dI <- n*γ*Iprev - (n*γ + μ)*I
    dR <- n*γ*I[[n]] - μ*R
    return(list(c(dS, dI, dR)))
  })
}

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
params <- c(β=0.1, γ = 1/12, μ = 0.001, n = n)
states <- c(0.9, 0.1, numeric(n-1), 0)
names(states) <- c("S", paste0("I", 1:n), "R") 

soln <- ode(y = states,
            times = time,
            func = SInR,
            parms = params)
plot(soln)
head(rowSums(soln[,-1]))

######################################################################

SEmInR <- function(time, states, params){
  with(as.list(c(params)), {
    S <- states[[1]]
    E <- states[2:(m+1)]
    I <- states[(m+2):(m+n+1)]
    R <- states[[m+n+2]]
    
    dS <- μ - β*S*sum(I) - μ*S
    Eprev <- c(β*S*sum(I), E[1:(m-1)])
    dE <- m*σ*Eprev - (m*σ + μ)*E
      
    Iprev <- c(σ*E[[m]], I[1:(n-1)])
    dI <- n*γ*Iprev - (n*γ + μ)*I
    
    dR <- n*γ*I[[n]] - μ*R
    return(list(c(dS, dE, dI, dR)))
  })
}

m = 3
n = 2
params <- c(β=0.1, γ = 1/12, σ = 0.1 ,μ = 0.001, n = n, m = m)
states <- c(0.9, numeric(m), 0.1, numeric(n-1), 0)
names(states) <- c("S", paste0("E", 1:m), paste0("I", 1:n), "R") 

soln <- ode(y = states,
            times = time,
            func = SEmInR,
            parms = params)
plot(soln)
head(rowSums(soln[,-1]))

