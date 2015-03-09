library(knitr)
library(markdown)

setwd("/home/Shared/data/array/Microarray_Edwin")

inputK <- "/home/gosia/R/R_Microarrays_Edwin/Analysis_Mouse.Rmd"

# knitr::knit2html(inputK)

knitr::knit(inputK, tangle = TRUE)

inputM <- "/home/Shared/data/array/Microarray_Edwin/Analysis_Mouse.md"

markdown::markdownToHTML(inputM)


# library(rmarkdown)
# rmarkdown::render("/home/gosia/R/R_Microarrays_Caroline/Analysis.Rmd")


