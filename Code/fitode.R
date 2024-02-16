library(fitode)
library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)

#######################################################################

# generate SInR formula (list)
sinrFG <- function(n) {
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
  
  fs <- purrr::map2(resp, vars, ~reformulate(.x, response = .y))
  fs <- lapply(fs, function(f) { environment(f) <- NULL; f })
  
  return(fs)
}

######################################################################

# generate initial value list
iniVG <- function(n) {
  Ivec <- paste0("I",1:n)
  vars <- c("S", Ivec, "R", "Cinc")
  inip <- c("N - i0",
            "i0",
            rep(0, n+1))
  
  Ini <- purrr::map2(inip, vars, ~reformulate(.x, response = .y))
  Ini <- lapply(Ini, function(f) { environment(f) <- NULL; f })
  
  return(Ini)
}

######################################################################

n <- 4

# include into a function ##
SIR_model <- odemodel(
  name="SInR (nbinom)",
  model=sinrFG(n),
  observation=list(
    confirmed ~ dnbinom(mu=Cinc, size=nbs)
  ),
  initial=iniVG(n),
  #diffnames="Cinc",
  par=c("beta", "gamma", "μ", "N", "i0", "nbs", "n"),
  #link=c(i0="logit")
)

######################################################################

β <- 0.2
D <- 10
μ <- 0.01  
N <- 10000
I0 <- 10
ts <- 1
T <- 200

nbs <- 1000
arp <- 0.9

######################################################################

SIR_start <- c(beta=β, gamma=1/D, μ = μ, N=10000, i0=I0, nbs=nbs, n = n)
ss_SIR <- simulate(SIR_model,
                   parms=SIR_start, times=timeSeq(1,200, FALSE))

######################################################################

sinr <- SInRFlow(β, D, n, μ, N, I0, ts, T)

df <- data.frame(Time = timeSeq(ts, T),
                 currentInc = diff(sinr[,'Cinc']),
                 fitodeInc = diff(ss_SIR[,'Cinc']))

######################################################################

print(ggplot(df, aes(x=Time))
      + geom_line(aes(y=currentInc, color="original inc"), linewidth = 1.5)
      + geom_line(aes(y=fitodeInc, color="fitode inc"), linetype = "dashed", linewidth = 2) 
      + labs(y = "Incidence", title = "Original Inc vs. FitOde Inc")
)



# change S0 to N
# set nbs in function instead of 1000
# change C inc in origin function


