library(knitr)
library(markdown)

setwd("/home/Shared/data/array/Microarray_Edwin")

inputK <- "/home/gosia/R/R_Microarrays_Edwin/Analysis_Mouse.Rmd"

knitr::knit(inputK, tangle = TRUE)
knitr::knit2html(inputK)



# library(rmarkdown)
# rmarkdown::render("/home/gosia/R/R_Microarrays_Edwin/Analysis.Rmd")





##########################################################################


library(knitr)
library(markdown)

setwd("/home/Shared/data/array/Microarray_Edwin")

inputK <- "/home/gosia/R/R_Microarrays_Edwin/Getting_probe_set_annotation.Rmd"

knitr::knit(inputK, tangle = TRUE)
knitr::knit2html(inputK)








