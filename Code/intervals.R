## Rounded intervals mean that we just choose numbers and round them

set.seed(223)

numObs <- 100
kappa <- 0.4
D <- 5

## k = 1/shape
## D = shape*scale = scale/kappa

v <- round(rgamma(numObs, shape=1/kappa, scale=D*kappa))

print(mean(v))
print(var(v)/mean(v)^2)

## start/end times

## start is a uniform
## end is start + v [above]
## round them both before subtracting
