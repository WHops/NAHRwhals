# Testscript
# 
# library(optparse)
# 
# 
# option_list <- list ( make_option (c("-f","--filelist"),default="blah.txt", 
#                                    help="comma separated list of files (default %default)")
# )
# 
# opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
# 
# myfilelist <- strsplit(opt$filelist, ",")
# 
# print(opt$filelist)
# print(myfilelist[[1]][1])