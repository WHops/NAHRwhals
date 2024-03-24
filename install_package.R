#!/usr/bin/env Rscript

# Re-install reshape2, which is messed up by conda currently (to be fixed properly asap)
#install.packages("reshape2", repos="https://cran.us.r-project.org")

# Use devtools to install the package. 
devtools::document()
#devtools::install()

