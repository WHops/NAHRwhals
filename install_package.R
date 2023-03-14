#!/usr/bin/env Rscript

# Install devtools if not installed.
if(!require(devtools)){
    install.packages("devtools")
}

# Use devtools to install the package. 
devtools::install()
