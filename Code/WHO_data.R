
library(readr)

library(shellpipes)

url <- "https://raw.githubusercontent.com/Tiffany-Xie/Thesis_NX/main/data/WHO-COVID-19-global-data.csv?token=GHSAT0AAAAAACN7IGL2VRCX4IX36MUIQKMMZQZNO2Q"

data <- read_csv(url)

saveVars(data)
