library(deSolve)
library(devtools)
#install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)

######################################################################

simObs <- function(sinr, arp, nbs, seed) {
  set.seed(seed)
  inc <- diff(sinr[,"inc"]) #/ts
  obs <- rnbinom(mu=arp*inc, size=nbs, n=length(inc))
  df <- data.frame(Time = (sinr[,"time"][-1] + sinr[,"time"][-dim(sinr)[1]])/2,
                   inc = inc,
                   obs = obs)
  return(df)
}

fitting <- function(df, βe, De, n, μ, S0, I0, ts, T, seed, plot=TRUE, optMethod="Nelder-Mead") {
  obs = df$obs
  trace_betae <<- c()
  trace_De <<- c()
  sir.nll <- function(βe, De, obs) {
    trace_betae <<- c(trace_betae, βe)
    trace_De <<- c(trace_De, De)
    
    out <- as.data.frame(SInRFlow(β=exp(βe), D=exp(De), n=n, μ=μ, S0=S0, I0=I0, ts=ts, T=T))
    nll <- -sum(dnbinom(x=obs, mu=diff(out$inc), size=1000, log=TRUE))
  }
  
  params0 <-list(βe=βe, De=De)
  set.seed(seed) 
  fit0 <- mle2(sir.nll, start=params0, data=list(obs=obs), method=optMethod,
               control = list(trace = 0))
  
  if (plot) {
    plotting(df, fit0, trace_betae, trace_De)
  }
  return(fit0)
  
}


plot_sim <- function(df, fit) {
  mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(coef(fit)[["βe"]]), D=exp(coef(fit)[["De"]]), n, μ, S0, I0, ts, T)))
  df["fitInc"] = diff(mod.prep$inc)
  
  p <- ggplot(df, aes(x=Time)) +
        geom_line(aes(y=obs, color='Observed')) +
        geom_line(aes(y=inc, color = 'Incidence'), linewidth=1) +
        geom_line(aes(y=fitInc, color = 'Fit Incidence'), linewidth=1.5, alpha=0.8)
  return(p)
}

plot_trace <- function(trace_betae, trace_De) {
  pardf <- data.frame(order = seq(0:(length(trace_betae)-1)), 
                      betae = trace_betae,
                      De = trace_De)
  
  trace1 <- ggplot(pardf, aes(x = De, y = betae)) +
    geom_point(aes(color = order)) +  # Color by Time
    geom_path(alpha = 0.5) +  # Trace path
    scale_color_gradient(low = "blue", high = "red") +  # Color gradient
    labs(x = "log(D)", y = "log(beta)", color = "Time")
  
  trace2 <- ggplot(pardf) +
        aes(x = order, y = De) +
        geom_point()
  
  trace3 <- ggplot(pardf) +
        aes(x = order, y = betae) +
        geom_point()
  
  return(list(trace1=trace1, trace2=trace2, trace3=trace3))
}

plotting <- function(df, fit0, trace_betae, trace_De) {
  p1 <- plot_sim(df, fit0)
  plots <- plot_trace(trace_betae, trace_De)
  p2 <- plots$trace1
  p3 <- plots$trace2
  p4 <- plots$trace3
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

#plotting <- function(df, fit, ts, T, trace_betae, trace_De) {
#  plot(timeSeq(ts, T), df$obs, type="l", xlab='Time, (Days)', ylab='I(t)')
#  t <- timeSeq(ts, T)
#  mod.prep <- as.data.frame(as.data.frame(SInRFlow(β=exp(coef(fit)[["βe"]]), D=exp(coef(fit)[["De"]]), n, μ, S0, I0, ts, T)))
#  lines(diff(mod.prep$inc)~t, col = "red", lwd=3)
#  lines(x=timeSeq(ts, T), y=df$inc, col = "blue")
#}

## Always works ####################################################################     

β <- 0.2
D <- 10
n <- 4
μ <- 0.01
S0 <- 999
I0 <- 1
ts <- 1
T <- 200

arp <- 0.9
nbs <- 1000

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, seed=72)
trace_betae <- c()
trace_De <- c()
ans <- fitting(df, βe=-1, De=1, n, μ, S0, I0, ts, T, seed=77
               , optMethod="Nelder-Mead"
)
print(ans)

pardf <- data.frame(order = seq(0:(length(trace_betae)-1)), 
                    betae = trace_betae,
                    De = trace_De)

ggplot(pardf, aes(x = De, y = betae)) +
  geom_point(aes(color = order)) +  # Color by Time
  geom_path(alpha = 0.5) +  # Trace path
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient
  labs(x = "D", y = "Beta", color = "Time", title = "success")

print(ggplot(pardf)
      + aes(x = order, y = De)
      + geom_point() # Color by Time
)

print(ggplot(pardf)
      + aes(x = order, y = betae)
      + geom_point() # Color by Time
)

## Change time #################################################################### Change time

T <- 200

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 71) # seed = 71

trace_betae <- c()
trace_De <- c()
tryCatch({
  fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, 72) # seed = 72
  1
}, warning = function(w) {
  pardf <- data.frame(order = seq(0:(length(trace_betae)-1)), 
                      betae = trace_betae,
                      De = trace_De)
  
  ggplot(pardf, aes(x = De, y = betae)) +
    geom_point(aes(color = order)) +  # Color by Time
    geom_path(alpha = 0.5) +  # Trace path
    scale_color_gradient(low = "blue", high = "red") +  # Color gradient
    labs(x = "D", y = "Beta", color = "Time", title = "failure")
}, error = function(e) {
  0
}, finally = {
  1
})



## Change population (N) #################################################################### Change population (N)

T <- 100
S0 <- 99999

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 70) # seed = 70
fit = fitting(df, βe=-1.5, De=2, n, μ, S0, I0, ts, T, 69) # seed = 69

## Change ts ####################################################################

S0 <- 999
ts <- 0.1 # works

sinr <- SInRFlow(β, D, n, μ, S0, I0, ts, T)
df = simObs(sinr, arp, nbs, 67) # seed = 70
fit = fitting(df, βe=-0.5, De=2, n, μ, S0, I0, ts, T, 66) # seed = 69


quit()
######################################################################









