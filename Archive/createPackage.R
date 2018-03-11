# install packages
# install.packages("devtools")
library("devtools")
# devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd("/Users/Hause/Dropbox/Working Datasets/functionsR/")
create("hauseFunctions")

setwd("./hauseFunctions")
document()
