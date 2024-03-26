library(ggplot2)
library(shellpipes)
startGraphics(height=1, width=1)


df <- data.frame(A = seq(1,10), B = seq(1,10)^2)

ggplot(df, aes(x = A, y = B)) +
geom_line() + 
labs(title = "THIS is a test plot (TESE)")

# not working qswl....
