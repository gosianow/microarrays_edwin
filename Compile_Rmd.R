library(knitr)
library(markdown)

setwd("/home/Shared/data/array/Microarray_Edwin")

inputK <- "/home/gosia/R/R_Microarrays_Edwin/Analysis_Mouse.Rmd"

# knitr::knit2html(inputK)

knitr::knit(inputK, tangle = TRUE)
knitr::knit(inputK, tangle = FALSE)

inputM <- "/home/Shared/data/array/Microarray_Edwin/Analysis_Mouse.md"

markdown::markdownToHTML(inputM)


# library(rmarkdown)
# rmarkdown::render("/home/gosia/R/R_Microarrays_Caroline/Analysis.Rmd")





##########################################################################


library(knitr)
library(markdown)

setwd("/home/Shared/data/array/Microarray_Edwin")

inputK <- "/home/gosia/R/R_Microarrays_Edwin/Getting_probe_set_annotation.Rmd"

knitr::knit(inputK, tangle = TRUE)
knitr::knit2html(inputK)








