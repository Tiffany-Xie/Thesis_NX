r2kappa <- function(r, n, offset=0){
	delta <- (1:n) - (n+1)/2
	res <- exp(delta*log(r))
	kappa <- sum(res^2)/(sum(res)^2)
	return(kappa - offset)
}

r2kappaF <- function(r, n, offset=0){
	return((r^n+1)*(r-1) / ((r^n-1)*(r+1)) - offset)
}

kappa2r <- function(kappa, n){
	rmax <- 7 ## BAD CODE, XNR please fix
	if(kappa>=1) return(NA)
	if(kappa<1/n) return(NA)
	u <- uniroot(r2kappa, interval=c(1, rmax), n=n, offset=kappa)
	return(u$root)
}

######################################################################

## Tests

r2kappa(1, 5)
r2kappa(2, 10)
r2kappaF(1, 5)
r2kappaF(2, 10)

kappa2r(1, 10)
kappa2r(0.1, 10)
kappa2r(0.4, 2)
