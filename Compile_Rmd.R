library(knitr)
library(rmarkdown)
setwd("/home/Shared/data/array/Microarray_Edwin")



inputK <- "/home/gosia/R/R_Microarrays_Edwin/Analysis_Mouse.Rmd"

knitr::knit(inputK, tangle = TRUE)
# knitr::knit2html(inputK)

rmarkdown::render(inputK, output_dir = "/home/Shared/data/array/Microarray_Edwin")





inputF <- "/home/gosia/R/R_Microarrays_Edwin/Plots_Figure1.Rmd"

rmarkdown::render(inputF, output_dir = "/home/Shared/data/array/Microarray_Edwin")





##########################################################################


# library(knitr)
# library(markdown)
# 
# setwd("/home/Shared/data/array/Microarray_Edwin")
# 
# inputK <- "/home/gosia/R/R_Microarrays_Edwin/Getting_probe_set_annotation.Rmd"
# 
# knitr::knit(inputK, tangle = TRUE)
# knitr::knit2html(inputK)








