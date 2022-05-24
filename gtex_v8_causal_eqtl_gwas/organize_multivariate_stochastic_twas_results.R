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
tissue_names <- read.table(gtex_tissue_file, header=TRUE, sep="\t")$pseudotissue_name



component_arr <- c()
tissue_arr <- c()
gene_arr <- c()
nominal_twas_z_score_arr <- c()
nominal_twas_pvalue_arr <- c()
nominal_stochastic_z_score_mean_arr <- c()
nominal_stochastic_z_score_sd_arr <- c()
multivariate_twas_pip_arr <- c()
multivariate_twas_stochastic_avg_pip_arr <- c()



for (chrom_num in 1:22) {
if (chrom_num != 6) {
multivariate_twas_organization_file <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_", chrom_num, "_component_stochastic_multivariate_twas_overlaps.txt")
multivariate_twas_df <- read.table(multivariate_twas_organization_file, header=TRUE)
num_components <- dim(multivariate_twas_df)[1]

for (component_num in 1:num_components) {
	# Extract relevent fields for component
	component_lead_variant_name <- multivariate_twas_df$lead_variant_name[component_num]
	component_window_id <- multivariate_twas_df$window_id[component_num]
	component_name <- paste0(component_lead_variant_name,"_", component_window_id)

	# TWAS data relevent to this component
	component_gene_tissue_names_file <- multivariate_twas_df$gene_tissue_names_file[component_num]
	component_nominal_z_score_obj <- multivariate_twas_df$nominal_z_score_obj[component_num]
	stochastic_nominal_z_score_obj <- multivariate_twas_df$stochastic_nominal_z_score_obj[component_num]
	component_pip_obj <- multivariate_twas_df$pip_obj[component_num]
	stochastic_pip_obj <- multivariate_twas_df$stochastic_pip_obj[component_num]

	if (!is.na(component_gene_tissue_names_file) && !is.na(component_nominal_z_score_obj) && !is.na(component_pip_obj) && !is.na(stochastic_pip_obj)) {

		# Contains column called ordered "gene" and "tissue info"
		gene_tissue_names <- readRDS(component_gene_tissue_names_file)
		tissue_names_arr <- c()
		gene_names_arr <- c()
		for (itera in 1:length(gene_tissue_names)) {
			gene_tissue_name <- gene_tissue_names[itera]
			info = strsplit(gene_tissue_name, ":")
			tissue_names_arr = c(tissue_names_arr, info[[1]][2])
			gene_names_arr = c(gene_names_arr, info[[1]][1])
		}

		# Extract nominal z-scores (and then extract pvalues)
		nom_z_scores <- readRDS(component_nominal_z_score_obj)
		nom_pvalues <- 2*(pnorm( abs(nom_z_scores) , lower.tail=F))

		# Extract stochastic z-scores
		stochastic_nominal_zscores <- readRDS(stochastic_nominal_z_score_obj)
		stochastic_zscore_mat <- matrix(unlist(stochastic_nominal_zscores), nrow = 100, byrow = TRUE)

		stochastically_averaged_zscores = colMeans(stochastic_zscore_mat, na.rm=TRUE)
		stochastically_estimated_zscore_standard_errors = apply(stochastic_zscore_mat,2,sd)

		# Extract PIPs
		pips <- readRDS(component_pip_obj)

		# Extract stochastic pips
		stochastic_pips <- readRDS(stochastic_pip_obj)
		stochastic_pips_mat <- matrix(unlist(stochastic_pips), nrow = 100, byrow = TRUE)
		stochastically_averaged_pips = colMeans(stochastic_pips_mat, na.rm=TRUE)


		# Add tp arrays to keep track
		component_arr <- c(component_arr, rep(component_name, length(tissue_names_arr)))
		tissue_arr <- c(tissue_arr, tissue_names_arr)
		gene_arr <- c(gene_arr, gene_names_arr)
		nominal_twas_z_score_arr <- c(nominal_twas_z_score_arr, nom_z_scores)
		nominal_stochastic_z_score_mean_arr <- c(nominal_stochastic_z_score_mean_arr, stochastically_averaged_zscores)
		nominal_stochastic_z_score_sd_arr <- c(nominal_stochastic_z_score_sd_arr, stochastically_estimated_zscore_standard_errors)
		nominal_twas_pvalue_arr <- c(nominal_twas_pvalue_arr, nom_pvalues)
		multivariate_twas_pip_arr <- c(multivariate_twas_pip_arr, pips)
		multivariate_twas_stochastic_avg_pip_arr <- c(multivariate_twas_stochastic_avg_pip_arr, stochastically_averaged_pips)
	}
}
}
}

df <- data.frame(trait_component=component_arr, tissue=tissue_arr, gene=gene_arr, nominal_twas_zscore=nominal_twas_z_score_arr, nominal_stochastic_zscore_mean=nominal_stochastic_z_score_mean_arr,nominal_stochastic_zscore_sd=nominal_stochastic_z_score_sd_arr,  nominal_twas_pvalue=nominal_twas_pvalue_arr, multivariate_twas_pip=multivariate_twas_pip_arr, multivariate_twas_stochastically_averaged_pip=multivariate_twas_stochastic_avg_pip_arr)

output_file <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

write.table(df,output_file, quote=FALSE, sep="\t",row.names = FALSE)