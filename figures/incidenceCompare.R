library(deSolve)
library(ggplot2); theme_set(theme_minimal(base_size = 12))
library(bbmle)
library(pseudoErlang)

library(shellpipes)
startGraphics(height=5, width=7)

######################################################################

β <- 0.2
D <- 10
fixn <- 12
μ <- 0.001
ts <- 0.1
T <- 100

N = 1
I0 = 0.01

n <- 4
kappa <- 1/n
time <- timeSeq(ts, T)
sinner <- sinnerFlow(β, D, kappa, fixn, μ, N=N, I0=I0, ts, T)
sinr <- SInRFlow(β, D, n, μ, N=N, I0=I0, ts, T)
dfn4 <- data.frame(Time = time, 
                   LCT = diff(sinr[,"Cinc"])/ts,
                   GCT = diff(sinner[,"Cinc"])/ts)
n <- 8
kappa <- 1/n
time <- timeSeq(ts, T)
sinner <- sinnerFlow(β, D, kappa, fixn, μ, N=N, I0=I0, ts, T)
sinr <- SInRFlow(β, D, n, μ, N=N, I0=I0, ts, T)
dfn8 <- data.frame(Time = time, 
                   LCT = diff(sinr[,"Cinc"])/ts,
                   GCT = diff(sinner[,"Cinc"])/ts)


df <- data.frame(Time = time, LCT4 = dfn4$LCT, LCT8 = dfn8$LCT,
                 GCT4 = dfn4$GCT, GCT8 = dfn8$GCT)

ggplot(df, aes(x=Time)) +
  geom_line(aes(y=LCT4, color = "n = 4"), linewidth=1) +
  geom_line(aes(y=LCT8, color = "n = 8"), linewidth=1) +
  scale_color_manual(values = c("n = 4" = "blue", "n = 8" = "orange")) +
  labs(color = "LCT", x = "Time, days",
       y = "Newly Infected Cases (proportion)") +
  geom_line(aes(y=GCT4, color = "black"), linewidth=1, linetype="dashed") +
  geom_line(aes(y=GCT8, color = "black"), linewidth=1, linetype="dashed")



