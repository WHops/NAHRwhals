#!/usr/bin/env Rscript

# Install devtools if not installed.
if(!require(devtools)){
    install.packages("devtools")
}

# Also this we will need for the nahrwhals.R wrapper
if(!require(argparse)){
    install.packages("argparse")
}


# Use devtools to install the package. 
devtools::install()
