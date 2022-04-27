args = commandArgs(trailingOnly=TRUE)
library(hash)
library(susieR)
options(warn=1)






run_causal_effect_size_regression <- function(trait_susie_file, tissue_susie_files, tissues) {
	# Extract predicted effect sizes from trait
	trait_res <- readRDS(trait_susie_file)
	beta_trait = coef(trait_res)
	beta_trait = as.numeric(beta_trait[2:length(beta_trait)]) # remove intercept

	num_trait_components = length(trait_res$sets$cs)


	# Number of snps in this analysis
	P = length(beta_trait)

	# Number of tissues
	T = length(tissue_susie_files)

	# Initilize tissue beta
	beta_tissue_mat <- matrix(0,nrow=P,ncol=T)
	tissue_components <- c()

	nominal_pvalues <- c()
	for (tt in 1:T) {
		tissue_name <- tissues[tt]
		tissue_susie_file <- tissue_susie_files[tt]
		tissue_res <- readRDS(tissue_susie_file)
		beta_tissue = coef(tissue_res)
		beta_tissue = as.numeric(beta_tissue[2:length(beta_tissue)]) # remove intercept
		beta_tissue_mat[,tt] = beta_tissue
		tissue_comp <- length(tissue_res$sets$cs)
		tissue_components <- c(tissue_components, tissue_comp) 
		lm_vers = lm(beta_trait~beta_tissue)
		coefs <- data.frame(coef(summary(lm_vers)))
		nominal_pvalues <- c(nominal_pvalues, coefs[2,4])
	}



	multivariate_pvalues <- rep(NA, T)
	pips <- rep(NA, T)

	#res <- susie(beta_tissue_mat, beta_trait, L=10, coverage=.99)
	#print(summary(res))
	if (num_trait_components > 0) {
		if (sum(tissue_components) > 0) {
			valid_tissues <- tissue_components > 0
			lm_vers = lm(beta_trait~beta_tissue_mat[,valid_tissues])
			coefs <- data.frame(coef(summary(lm_vers)))
			pvalues <- coefs[,4]
			pvalues <- pvalues[2:length(pvalues)]
			multivariate_pvalues[valid_tissues] = pvalues
			#print(summary(lm_vers))

			if (sum(valid_tissues) > 1) {
				res <- susie(beta_tissue_mat[,valid_tissues], beta_trait, L=10)
				#print(summary(res))
				pips[valid_tissues] <- res[['pip']]
				#print(res[['pip']])
			}

		}
	}
	return(list(num_trait_components=num_trait_components, num_tissue_components=tissue_components, nominal_pvalues=nominal_pvalues, multivariate_pvalues=multivariate_pvalues, pips=pips))
}







trait_name = args[1]
gene_file = args[2]
gtex_susie_results_dir = args[3]
ukbb_susie_results_dir = args[4]
gtex_pseudotissue_file = args[5]
output_dir = args[6]
job_number = as.numeric(args[7])
total_jobs = as.numeric(args[8])




full_gtex_tissue_df <- read.table(gtex_pseudotissue_file, header=TRUE, sep="\t")
full_tissue_names <- full_gtex_tissue_df$pseudotissue_name

num_tissues = length(full_tissue_names)

tissue_name_to_position <- hash()
for (tissue_iter in 1:length(full_tissue_names)) {
	tissue_name_to_position[[full_tissue_names[tissue_iter]]] = tissue_iter
}



# Load in gene data frame
gene_df <- read.table(gene_file, header=TRUE, sep="\t")
# Get total number of genes
num_genes <- dim(gene_df)[1]


# For parallelization purposes, determine which genes to test in this thread
tasks_per_job = floor(num_genes/total_jobs) + 1
start_task = floor(job_number*tasks_per_job + 1)
end_task = floor((job_number + 1)*tasks_per_job)
if (end_task > num_genes) {
	end_task = num_genes
}


output_file <- paste0(output_dir, "causal_effect_regression_", job_number, "_", total_jobs, ".txt")
sink(output_file)

liner <- "gene_name\ttissue\tnum_trait_components\tnum_tissue_components\tnominal_pvalue\tmultivariate_pvalue\tpip\n"
cat(liner)


for (gene_num in start_task:end_task) {
	gene_name <- gene_df$gene_id[gene_num]
	variant_file <- gene_df$variant_file[gene_num]
	tissue_file <- gene_df$tissue_file[gene_num]

	variant_ids <- read.table(variant_file)$V1
	tissues <- read.table(tissue_file)$V1

	# Trait Susie file
	trait_susie_file <- paste0(ukbb_susie_results_dir, gene_name, "_", trait_name, "_susie_res.RDS")

	tissue_susie_files <- c()
	used_tissues <- c()
	for (tissue_iter in 1:length(tissues)) {
		tissue_name <- tissues[tissue_iter]
		gtex_tissue_susie_file <- paste0(gtex_susie_results_dir, gene_name, "_", tissue_name, "_susie_res.RDS")
		if (file.exists(gtex_tissue_susie_file)) {
			tissue_susie_files <- c(tissue_susie_files, gtex_tissue_susie_file)
			used_tissues <- c(used_tissues, tissues[tissue_iter])
		}
	}

	if (length(tissue_susie_files) > 0) {
		if (file.exists(trait_susie_file)) {
			# At this point we can confirm susie converged for this gene
			causal_effects_res <- run_causal_effect_size_regression(trait_susie_file, tissue_susie_files, used_tissues)
			for (itera in 1:length(used_tissues)) {
				liner <- paste0(gene_name, "\t", used_tissues[itera], "\t", causal_effects_res[['num_trait_components']], "\t", causal_effects_res[['num_tissue_components']][itera], "\t", causal_effects_res[['nominal_pvalues']][itera], "\t",causal_effects_res[['multivariate_pvalues']][itera], "\t",causal_effects_res[['pips']][itera] , "\n")
				cat(liner)
			}

		}
	}


}

sink()
