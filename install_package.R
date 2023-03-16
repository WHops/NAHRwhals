#!/usr/bin/env Rscript

# Install devtools if not installed.
if(!require(devtools)){
    install.packages("devtools", repos="https://cran.us.r-project.org")
}

# Also this we will need for the nahrwhals.R wrapper
if(!require(argparse)){
    install.packages("argparse", repos="https://cran.us.r-project.org")
}

# Use devtools to install the package. 
devtools::install()
