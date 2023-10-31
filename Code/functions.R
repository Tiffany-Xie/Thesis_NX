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
	u <- uniroot(r2kappa, interval=c(1, 10), n=n, offset=kappa)
	return(u$root)
}

######################################################################

## Tests

r2kappa(1, 5)
r2kappa(2, 10)
r2kappaF(1, 5)
r2kappaF(2, 10)

k <- kappa2r(0.4, 10)
print(k)
r2kappa(k, 10)
