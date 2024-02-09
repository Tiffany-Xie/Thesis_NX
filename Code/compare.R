
library(deSolve)
library(ggplot2); theme_set(theme_minimal())
library(bbmle)

## install_github("Tiffany-Xie/pseudoErlang")
library(pseudoErlang)

β <- 0.2
D <- 10
n <- 4
μ <- 0.01  
S0 <- 999
I0 <- 1
ts <- 1
T <- 200

sinr1 <- SInRFlow(β, D, n=1, μ, S0, I0, ts, T)
sinr2 <- SInRFlow(β, D, n=2, μ, S0, I0, ts, T)
sinr4 <- SInRFlow(β, D, n=4, μ, S0, I0, ts, T)
sinr8 <- SInRFlow(β, D, n=8, μ, S0, I0, ts, T)
time <- timeSeq(ts, T)
plot(time, diff(sinr8[,"inc"]), 
     main = "Erlang comparison", 
     type="l")
lines(time, diff(sinr4[,"inc"]), col="orange")
lines(time, diff(sinr2[,"inc"]), col="blue")
## lines(time, diff(sinr1[,"inc"]), col="red")
