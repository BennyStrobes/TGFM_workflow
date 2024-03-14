args = commandArgs(trailingOnly=TRUE)
library(ctwas)
options(warn=1)







#####################
# Command line args
#####################
chrom_num = args[1]
processed_ctwas_genotype_data_dir = args[2]
ctwas_regions_file = args[3]

# Load in quasi-independent regions from causal-twas data
regions <- system.file("extdata/ldetect", "EUR.b37.bed", package = "ctwas")
regions_df <- read.table(regions, header = T)

# Filter to only regions on chromosome 1
regions_df_subset <- regions_df[as.character(regions_df$chr) == paste0("chr", chrom_num),]

# Print to output
write.table(regions_df_subset, file=ctwas_regions_file, row.names=F, col.names=T, sep="\t", quote = F)

