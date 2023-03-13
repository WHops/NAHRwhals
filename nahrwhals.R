#!/bin/Rscript

# Load nahrwhals and 
library(devtools)
devtools::load_all()
#library(nahrtoolkit)
library(argparse)

###### INPUT: CONFIG plus CMDLINE-OVERWRITING ########

# Define command line arguments
parser <- ArgumentParser()
parser$add_argument("--config-file", default="config.txt", help="Configuration file path")
parser$add_argument("--params", nargs="*", help="Override configuration parameters")
args <- parser$parse_args()

# Load configuration file
config <- read.table(args$config_file, sep="=", strip.white=TRUE, comment.char = "#")
config <- setNames(config$V2, config$V1)

# Override configuration file with command line arguments
for (arg in args$params) {
  arg_name <- strsplit(arg, "=")[[1]][1]
  arg_value <- strsplit(arg, "=")[[1]][2]
  if (!(arg_name %in% names(config))) {
    stop(paste0("Invalid argument '", arg_name, "'. '", arg_name, "' is not a valid parameter name."))
  }
  config[[arg_name]] <- arg_value
}

# Store configuration parameters in a list
params <- list()
for (param_name in names(config)) {
  params[[param_name]] <- (config[[param_name]])
}

# Print parameter values
for (param_name in names(params)) {
  cat(param_name, "=", params[[param_name]], "\n")
}


boolean_characters = c("T", "TRUE", "F", "FALSE")
# Characters to numeric, booleans or NULL
for (i in seq_along(params)) {
  if (!((is.na(as.numeric(params[[i]]))))) {
    params[[i]] <- as.numeric(params[[i]])
  } else if (params[[i]] %in% boolean_characters) {
    params[[i]] <- as.logical(params[[i]])
  } else if (is.null(params[[i]])) {
    params[[i]] <- NULL
  } 
}



###### Compute parameters that are functions of other parameters. #####
default_param_values = c('auto', 'Auto', '', NA, NULL, F)

if (params$chunklen %in% default_param_values){
  params$chunklen = determine_chunklen(params$start_x, params$end_x)
}
if (params$minlen %in% default_param_values){
  params$minlen = determine_compression(params$start_x, params$end_x)
}
if (params$compression %in% default_param_values){
  params$compression = determine_compression(params$start_x, params$end_x)
}

# This one honestly I can't remember what it is let's just keep it
if (params$baseline_log_minsize_max %in% default_param_values){
  params$baseline_log_minsize_max = max(log2(20000), log2((params$end_x-params$start_x) / 10))
}



##### Run the main NAHRwhals wrapper function #####
wrapper_aln_and_analyse(params)
  
print('done!')
