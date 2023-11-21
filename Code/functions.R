r2kappa <- function(r, n, offset=0){
	delta <- (1:n) - (n+1)/2
	res <- exp(delta*log(r))
	kappa <- sum(res^2)/(sum(res)^2)
	return(kappa - offset)
}

## Formula version probably not needed
r2kappaF <- function(r, n, offset=0){
	return((r^n+1)*(r-1) / ((r^n-1)*(r+1)) - offset)
}

kappa2r <- function(kappa, n){
	rmax <- 2*(1+kappa)/(1-kappa)
	if(kappa>=1) return(NA)
	if(kappa<1/n) return(NA)
	u <- uniroot(r2kappa, interval=c(1, rmax), n=n, offset=kappa)
	return(u$root)
}

######################################################################

gammaFlow <- function (mu, kappa, ts, T){
	steps <- ceiling(T/ts)
	boundaries <- seq(0, by=ts, length.out=steps+1)
	cum <- pgamma(boundaries, shape=1/kappa, scale=mu*kappa)
	return(diff(cum))
}

flowMoments <- function(flow, ts, offset=1/2){
	tot <- sum(flow)
	nflow <- flow/tot
	times <- ((1:length(flow))-offset)*ts
	mu <- sum(times*nflow)
	S <-  sum(times^2*nflow)
	kap <- S/mu^2-1
	return(c(tot=tot, mu=mu, kap=kap))
}

######################################################################

boxIntegrate <- function(time, gam, model) {
	n <- length(gam)
	names(states) <- c(paste0("I", 1:n), "R")
	cum <- ode(y = states,
		y=c(1, numeric(n))
		, times = time
		, func = model,
		parms = c(n = n, gam=gam)
	)
	return(cum)
}

######################################################################

quit()

ts <- 0.01
g <- gammaFlow(mu=5, kappa=0.4, ts=ts, T=50)
print(flowMoments(g, ts=ts, offset=1/2))

g <- gammaFlow(mu=5, kappa=1, ts=ts, T=50)
print(flowMoments(g, ts=ts, offset=1/2))

g <- gammaFlow(mu=5, kappa=1, ts=ts, T=100)
print(flowMoments(g, ts=ts, offset=1/2))

######################################################################

quit()

## Test kappa2r and relatives

r2kappa(1, 5)
r2kappa(2, 10)
r2kappaF(1, 5)
r2kappaF(2, 10)

kappa2r(0.2, 10)
r2kappa(kappa2r(0.2, 10), 10)
