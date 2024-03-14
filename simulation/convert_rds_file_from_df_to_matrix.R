args = commandArgs(trailingOnly=TRUE)
options(warn=1)













################
# Command line args
################
input_file <- args[1]
output_file <- args[2]



data <- as.matrix(readRDS(file=input_file))
saveRDS(data, file=output_file)
