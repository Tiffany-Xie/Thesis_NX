library(deSolve)
library(devtools)
#install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(ggplot2); theme_set(theme_bw(base_size=14))

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
