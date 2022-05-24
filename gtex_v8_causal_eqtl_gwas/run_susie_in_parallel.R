args = commandArgs(trailingOnly=TRUE)
library(hash)
library(susieR)
options(warn=1)






susie_run_across_traits_shell <- function(gene_name, variant_ids, trait_names, beta_mat, std_err_mat, LD, sample_sizes, output_dir) {
	for (trait_num in 1:length(trait_names)) {
		trait_name <- trait_names[trait_num]
		
		#res = susie_suff_stat(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]),R=LD, n=sample_sizes[trait_num])
		#res = susie_rss(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]),R=LD, n=sample_sizes[trait_num])
		tryCatch(
		{
			#res = susie_suff_stat(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]),R=LD, n=sample_sizes[trait_num])
			res = susie_rss(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]),R=LD, n=sample_sizes[trait_num])
			output_file <- paste0(output_dir, gene_name, "_", trait_name, "_no_ambiguous_variants_susie_res.RDS")
			saveRDS(res, file=output_file)
		 },
    		error = function(e) {
	    	print(paste0(gene_name, " ", trait_name, " errror"))
   		}
    	)
	}
}











########################
# Command line args
########################
gene_file = args[1]
output_dir = args[2]
job_number = as.numeric(args[3])
total_jobs = as.numeric(args[4])


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
print(start_task)
print(end_task)

for (gene_num in start_task:end_task) {
		print(gene_num)
		gene_name <- gene_df$gene_id[gene_num]
		genotype_file <- gene_df$ref_genotype_file[gene_num]
		beta_file <- gene_df$beta_file[gene_num]
		std_err_file <- gene_df$std_err_file[gene_num]
		variant_file <- gene_df$variant_file[gene_num]
		tissue_file <- gene_df$tissue_file[gene_num]
		sample_size_file <- gene_df$sample_size_file[gene_num]

		beta_mat <- read.table(beta_file)
		std_err_mat <- read.table(std_err_file)

		variant_ids <- read.table(variant_file)$V1
		tissues <- read.table(tissue_file)$V1
		sample_sizes <- read.table(sample_size_file)$V1

		geno_df <- read.table(genotype_file, header=FALSE, sep="\t")
		LD <- cor(geno_df)

		susie_run_across_traits_shell(gene_name, variant_ids, tissues, beta_mat, std_err_mat, LD, sample_sizes, output_dir)

}
