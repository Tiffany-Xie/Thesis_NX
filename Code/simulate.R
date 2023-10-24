# Packages ####
library(deSolve)
library(gridExtra)
library(dplyr)

library(shellpipes)

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

PEdens <- function(r, nPE, γ=0.1, maxTime=40, ts=0.1, μ=0)
{	
	states <- c(1, numeric(nPE))
	names(states) <- c(paste ("I", 1:nPE, sep=""), "R")
	time <- seq(0, maxTime, by=ts)

	params <- c(γ = γ, μ = μ, n = nPE, a = nPE*γ, r = r)
	soln <- ode(y = states,
		times = time, 
		func = SInR_geom, 
		parms = params
	)
	return(data.frame(
		t = (time[-1] + time[-length(time)])/2
		, d = diff(soln[, "R"]/ts)
	))
}

ts <- 0.1
df <- PEdens(2, 12, ts=ts)
print(df)

check <- sum(df$d*ts)
mean <- sum(df$t*df$d*ts)
S <- sum(df$t^2*df$d*ts)

print(c(check, mean, S))
print(S/mean^2-1)

