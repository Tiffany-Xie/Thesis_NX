library(deSolve)

# centered geometric series
cgs <- function(r, n){
	delta <- (1:n) - (n+1)/2
	return(exp(delta*log(r)))
}

r2kappa <- function(r, n, offset=0){
	res <- cgs(r, n)
	kappa <- sum(res^2)/(sum(res)^2)
	return(kappa - offset)
}

## Formula version probably not needed
r2kappaF <- function(r, n, offset=0){
	return((r^n+1)*(r-1) / ((r^n-1)*(r+1)) - offset)
}

kappa2r <- function(kappa, n){
	if(kappa>=1) return(NA)
	if(kappa<1/n) return(NA)
	rmax <- 2*(1+kappa)/(1-kappa)
	u <- uniroot(r2kappa, interval=c(1, rmax), n=n, offset=kappa)
	return(u$root)
}

######################################################################

boundaries <- function(ts, T){
	steps <- ceiling(T/ts)
	return(seq(0, by=ts, length.out=steps+1))
}

centers <- function(ts,T){
	b <- boundaries(ts, T)
	return((b[-1]+b[-length(b)])/2)
}

######################################################################

## XNR: it might be fun (but low priority) to write a gammaFlow function called gammaFlowDens, that calculates in a different way and compares precision

gammaFlow <- function (mu, kappa, ts, T){
	b <- boundaries(ts, T)
	cum <- pgamma(b, shape=1/kappa, scale=mu*kappa)
	return(diff(cum))
}

## This is not using boundaries/centers yet (it relies on length of flow, which seems OK)
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

cohortModel <- function(t, states, params) {
	with(params, {
		stopifnot(length(gam)==n)
		I <- states[1:n]
		R <- states[[n+1]]

		flow <- gam*I
		inflow <- c(0, flow[1:(n-1)])
		dI <- inflow-flow
		dR <- flow[[n]]
		return(list(c(dI, dR)))
	})
}

boxIntegrate <- function(time, gam, model=cohortModel) {
	n <- length(gam)
	states <- c(1, numeric(n))
	names(states) <- c(paste0("I", 1:n), "R")
	ser <- ode(y = states
		, times = time
		, func = model
		, parms = list(n = n, gam=gam)
	)
	cum <- as.data.frame(ser)[["R"]]
	return(diff(cum))
}

## Do an Erlang-like model and get the flow
flatFlow <- function(mu, n, ts, T, model=cohortModel){
	gamma = n/mu
	gam = rep(gamma, n)
	time = boundaries(ts, T)
	return(boxIntegrate(time, gam))
}

geomFlow <- function(mu, n, kappa, ts, T, model=cohortModel){
	r <- kappa2r(kappa, n)
	D = cgs(r, n)
	D = D*mu/(sum(D))
	time = boundaries(ts, T)
	return(boxIntegrate(time, gam=1/D))
}

######################################################################

ts <- 0.1
n <- 4
T <- 200
mu <- 5
pen <- 12

g <- gammaFlow(mu=mu, kappa=1/n, ts=ts, T=T)
print(flowMoments(g, ts=ts))

f <- flatFlow(mu=mu, n=n, ts=ts, T=T)
print(flowMoments(f, ts=ts))

m <- geomFlow(mu=mu, n=pen, kappa=1/n, ts=ts, T=T)
print(flowMoments(m, ts=ts))

######################################################################

c <- centers(ts, T)

plot(c, g, xlim=c(0, 25))
plot(c, f, xlim=c(0, 25))
plot(c, m, xlim=c(0, 25))

######################################################################

quit()

## Test kappa2r and relatives

r2kappa(1, 5)
r2kappa(2, 10)
r2kappaF(1, 5)
r2kappaF(2, 10)

kappa2r(0.2, 10)
r2kappa(kappa2r(0.2, 10), 10)

######################################################################

boundaries(0.2, 3)
centers(0.2, 3)

######################################################################
