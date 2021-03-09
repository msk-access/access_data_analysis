#!/usr/bin/env Rscript

library(knitr)
library(rmarkdown)

args=commandArgs(trailingOnly=TRUE)
render(args[1])
