# fitode
library(fitode)


n <- 4
Ivec <- paste0("I",1:n)
incidence <- sprintf("beta*S*(%s)/N", paste(Ivec, collapse = "+"))
Itrans <- function(i, death=FALSE) {
  if (death) return(sprintf("(n*gamma + μ)*I%d", i))
  return(sprintf("n*gamma*I%d", i))
}
vars <- c("S", Ivec, "R", "Cinc")
resp <- c(sprintf("μ*N - %s - μ*S", incidence),
          sprintf("%s - (n*gamma + μ)*I1", incidence),
          sprintf("%s - %s", Itrans((2:n)-1), Itrans(2:n, TRUE)),
          sprintf("%s - μ*R", Itrans(n)),
          incidence)


## could also use base-R Map(), or a for loop ...
fs <- purrr::map2(resp, vars, ~reformulate(.x, response = .y))
fs <- lapply(fs, function(f) { environment(f) <- NULL; f })
## initial value
inip <- c("N - i0",
          "i0",
          rep(0, n+1))

Ini <- purrr::map2(inip, vars, ~reformulate(.x, response = .y))
Ini <- lapply(Ini, function(f) { environment(f) <- NULL; f })

#####
SIR_model <- odemodel(
  name="SIR (nbinom)",
  model=fs,
  observation=list(
    confirmed ~ dnbinom(mu=Cinc, size=phi)
  ),
  initial=Ini,
  diffnames="Cinc",
  par=c("beta", "gamma", "μ", "N", "i0", "phi", "n"),
  #link=c(i0="logit")
)

#######

SIR_start <- c(beta=0.2, gamma=0.1, μ = 0.01, N=10000, i0=10, phi=1000, n = n)
ss_SIR <- simulate(SIR_model,
                   parms=SIR_start, times=timeSeq(1,200, FALSE))
head(ss_SIR)
plot(ss_SIR$times, ss_SIR$Cinc, type="l")
lines(df$Time, df$inc)

###################3
SIR_model <- odemodel(
  name="SIR (nbinom)",
  model=list(
    S ~ - beta * S * I/N,
    I ~ beta * S * I/N - gamma * I,
    R ~ gamma * I
  ),
  observation=list(
    confirmed ~ dnbinom(mu=R, size=phi)
  ),
  initial=list(
    S ~ N * (1 - i0),
    I ~ N * i0,
    R ~ 0
  ),
  diffnames="R",
  par=c("beta", "gamma", "N", "i0", "phi"),
  #link=c(i0="logit")
)
SIR_start <- c(beta=70, gamma=60, N=40000, i0=0.0004, phi=6)
SierraLeone2014b <- rbind(
  c(times=SierraLeone2014$times[1] -
      diff(SierraLeone2014$times)[1], confirmed=NA),
  SierraLeone2014
)
ss_SIR <- simulate(SIR_model,
                   parms=SIR_start, times=SierraLeone2014b$times)
plot(ss_SIR$times, ss_SIR$S, type="l")

plot(SierraLeone2014)
lines(ss_SIR$times, ss_SIR$R)


