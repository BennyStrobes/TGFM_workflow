args = commandArgs(trailingOnly=TRUE)
library(hash)
library(susieR)
options(warn=1)






susie_run_across_traits_shell <- function(gene_name, variant_ids, trait_names, beta_mat, std_err_mat, LD, sample_sizes, output_dir) {
	for (trait_num in 1:length(trait_names)) {
		tryCatch (
		{
			trait_name <- trait_names[trait_num]
			print(trait_name)
			#res = susie_suff_stat(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]),R=LD, n=sample_sizes[trait_num])
			res <- susie_rss(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]), R=LD, n=sample_sizes[trait_num])
			res[["variant_ids"]] = variant_ids
			output_file <- paste0(output_dir, gene_name, "_", trait_name, "_susie_res.RDS")
			saveRDS(res, file=output_file)
    	},
    	error = function(e) {
    		print("ERORORO")
       		print(window_num)
    	}
    	)
	}
}











########################
# Command line args
########################
window_file = args[1]
output_dir = args[2]
job_number = as.numeric(args[3])
total_jobs = as.numeric(args[4])

# Load in gene data frame
window_df <- read.table(window_file, header=TRUE, sep="\t")
# Brief tweak
window_df$window_id = paste0(window_df$chrom_num,":", window_df$start_pos_inclusive, ":", window_df$end_position_exclusive)
# Get total number of genes
num_windows <- dim(window_df)[1]

# For parallelization purposes, determine which genes to test in this thread
tasks_per_job = floor(num_windows/total_jobs) + 1
start_task = floor(job_number*tasks_per_job + 1)
end_task = floor((job_number + 1)*tasks_per_job)
if (end_task > num_windows) {
	end_task = num_windows
}


print(start_task)
print(end_task)
for (window_num in start_task:end_task) {
	print(window_num)


	window_id <- window_df$window_id[window_num]
	genotype_file <- window_df$ref_genotype_file[window_num]
	beta_file <- window_df$beta_file[window_num]
	std_err_file <- window_df$std_err_file[window_num]
	variant_file <- window_df$variant_file[window_num]
	tissue_file <- window_df$tissue_file[window_num]
	sample_size_file <- window_df$sample_size_file[window_num]

	beta_mat <- read.table(beta_file)
	std_err_mat <- read.table(std_err_file)

	variant_ids <- read.table(variant_file)$V1
	tissues <- read.table(tissue_file)$V1
	sample_sizes <- read.table(sample_size_file)$V1

	geno_df <- read.table(genotype_file, header=FALSE, sep="\t")
	LD <- cor(geno_df)

	susie_run_across_traits_shell(window_id, variant_ids, tissues, beta_mat, std_err_mat, LD, sample_sizes, output_dir)

}

