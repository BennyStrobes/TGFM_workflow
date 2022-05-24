args = commandArgs(trailingOnly=TRUE)
library(hash)
library(susieR)
library('plink2R')
library("susieR")
options(warn=1)








trait_name <- args[1]
gtex_fusion_multivariate_associations_dir <- args[2]
ukbb_genome_wide_susie_organized_results_dir <- args[3]
gtex_tissue_file <- args[4]

# Extract susie genome wide trait components on this chromosome
organized_susie_trait_component_file <- paste0(ukbb_genome_wide_susie_organized_results_dir, trait_name, "_organized_susie_components.txt")  # File name
organized_susie_trait_component_df <- read.table(organized_susie_trait_component_file, header=TRUE)  # load in file

# Extract tissue names
tissue_names <- read.table(gtex_tissue_file, header=TRUE, sep="\t")$tissue_name 




component_arr <- c()
tissue_arr <- c()
gene_arr <- c()
nominal_twas_pvalue_arr <- c()
nominal_twas_zscore_arr <- c()
multivariate_twas_pip_arr <- c()


chrom_num="21"

for (chrom_num in 1:22) {
multivariate_twas_organization_file <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_", chrom_num, "_component_multivariate_twas_overlaps.txt")
multivariate_twas_df <- read.table(multivariate_twas_organization_file, header=TRUE)
num_components <- dim(multivariate_twas_df)[1]

for (component_num in 1:num_components) {
	# Extract relevent fields for component
	component_lead_variant_name <- multivariate_twas_df$lead_variant_name[component_num]
	component_window_id <- multivariate_twas_df$window_id[component_num]
	component_name <- paste0(component_lead_variant_name,"_", component_window_id)

	# TWAS data relevent to this component
	component_twas_nominal_association_table_file <- multivariate_twas_df$twas_nominal_association_table[component_num]
	component_twas_multivariate_association_obj <- multivariate_twas_df$twas_multivariate_twas_obj[component_num]

	if (!is.na(component_twas_nominal_association_table_file)) {
	if (!is.na(component_twas_multivariate_association_obj)) {

	# Contains column called "gene", "tissue", "twas_nominal_z", and "twas_nominal_p"
	twas_nominal_df <- read.table(component_twas_nominal_association_table_file, header=TRUE)

	# Load in SuSiE multivariate twas object
	twas_multivariate_obj <- readRDS(component_twas_multivariate_association_obj)
	#print(component_twas_multivariate_association_obj)
	#print(summary(twas_multivariate_obj))
	#print(twas_multivariate_obj$pip)
	#print(twas_multivariate_obj$sets$cs_index)
	#print(twas_multivariate_obj$converged)
	#1.0 - prod(1-aa$alpha[aa$V > 0,1])
	#num_var = dim(twas_multivariate_obj$alpha)[2]
	#new_pips <- rep(0, num_var)
	#for (var_num in 1:num_var) {
	#	if (length(twas_multivariate_obj$sets$cs_index) > 0) {
	#		new_pip <- 1.0 - prod(1.0 - twas_multivariate_obj$alpha[twas_multivariate_obj$sets$cs_index, var_num])
	#		new_pips[var_num] = new_pip
	#	}
	#}

	component_arr <- c(component_arr, rep(component_name, length(twas_nominal_df$tissue)))
	tissue_arr <- c(tissue_arr, twas_nominal_df$tissue)
	gene_arr <- c(gene_arr, twas_nominal_df$gene)
	nominal_twas_pvalue_arr <- c(nominal_twas_pvalue_arr, twas_nominal_df$twas_nominal_p)
	nominal_twas_zscore_arr <- c(nominal_twas_zscore_arr, twas_nominal_df$twas_nominal_z)
	multivariate_twas_pip_arr <- c(multivariate_twas_pip_arr, twas_multivariate_obj$pip)
	#multivariate_twas_pip_arr <- c(multivariate_twas_pip_arr, new_pips)
	}
	}
}

}

df <- data.frame(trait_component=component_arr, tissue=tissue_arr, gene=gene_arr, nominal_twas_zscore=nominal_twas_zscore_arr, nominal_twas_pvalue=nominal_twas_pvalue_arr, multivariate_twas_pip=multivariate_twas_pip_arr)

output_file <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

write.table(df,output_file, quote=FALSE, sep="\t",row.names = FALSE)