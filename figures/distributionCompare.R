library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 17))
library(bbmle)
library(pseudoErlang)
library(gridExtra) #

library(shellpipes)
startGraphics(height=5, width=7)

######################################################################

compare <- function(D, n, nfix, ts, T) {
  time <- timeSeq(ts, T)
  sddLCT <- erlang(time, n, 1/D)
  sigr <- sinnerFlow(β=0, D=10, kappa=1/n, n=nfix, μ=0, N=1, I0=1, ts=0.1, T=60)
  df <- data.frame(Time = time,
                   ssdLCT = sddLCT,
                   ssdGCT = diff(sigr[,"R"])/ts)
  ggplot(df, aes(x=Time)) +
    geom_line(aes(y=ssdLCT, color="LCT"), linewidth=1) +
    geom_line(aes(y=ssdGCT, color="GCT"), linewidth=1) +
    labs(x="Stage duration, days", y="Probability Density",
         title = paste0("n = ", n))
}

######################################################################

ts <- 0.1
T <- 60
D <- 10
nfix <- 12
nn <- c(2,4,6,8)

plots <- lapply(nn, function(nn) compare(D, nn, nfix, ts, T))
grid.arrange(grobs = plots, ncol = 2, nrow = 2)
