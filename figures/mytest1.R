library(ggplot2)
library(shellpipes)
startGraphics()

df <- data.frame(A = seq(1,10), B = seq(1,10)^2)

ggplot(df, aes(x = A, y = B)) +
geom_line() + 
labs(title = "THIS is a test plot")

# not working qswl....
