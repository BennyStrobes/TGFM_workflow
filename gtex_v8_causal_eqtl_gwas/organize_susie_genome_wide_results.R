args = commandArgs(trailingOnly=TRUE)
library(hash)
library(susieR)
options(warn=1)











########################
# Command line args
########################
window_file = args[1]
input_dir = args[2]
output_dir = args[3]



# Load in gene data frame
window_df <- read.table(window_file, header=TRUE, sep="\t")
# Brief tweak
window_df$window_id = paste0(window_df$chrom_num,":", window_df$start_pos_inclusive, ":", window_df$end_position_exclusive)
# Get total number of genes
num_windows <- dim(window_df)[1]

# Get traits
trait_file <- window_df$tissue_file[1]
traits <- read.table(trait_file)$V1



for (trait_num in 1:length(traits)) {
	trait_name <- traits[trait_num]
	print(trait_name)

	trait_output_file <- paste0(output_dir, trait_name, "_organized_susie_components.txt")
	sink(trait_output_file)

	header <- paste0("chrom_num\tlead_variant_position\tlead_variant_name\tlead_variant_pi\twindow_id\tcs_variants\tcs_pi\n")
	cat(header)

	for (window_num in 1:num_windows) {

		window_id <- window_df$window_id[window_num]
		# Get chrom_num
		chrom_num <- window_df$chrom_num[window_num]
		# Extract whether window is at beginning or end of the chromosome
		chrom_first_bool = window_df$chrom_first_window_boolean[window_num]
		chrom_last_bool = window_df$chrom_last_window_boolean[window_num]
		# Get window start and end positions
		start_position_inclusive <- window_df$start_pos_inclusive[window_num]
		end_position_exclusive <- window_df$end_position_exclusive[window_num]
		# Get window specific start and end position
		window_specific_start_position = start_position_inclusive + 1000000
		window_specific_end_position = start_position_inclusive + 2000000
		if (chrom_first_bool == "True") {
			window_specific_start_position = window_specific_start_position - 1000000
		}
		if (chrom_last_bool == "True") {
			window_specific_end_position = window_specific_end_position + 1000000
		}
		susie_res_file <- paste0(input_dir, window_id, "_", trait_name, "_susie_res.RDS")
		susie_res <- readRDS(susie_res_file)
		# Numnber of credible sets passing our coverage threshold
		num_cs = length(susie_res$sets$cs)
		if (num_cs > 0) {

			for (cs_iter in 1:num_cs) {

				cs_num = susie_res$sets$cs_index[cs_iter]
				cs_variant_indices = susie_res$sets$cs[[cs_iter]]
				cs_variant_ids = susie_res$variant_ids[susie_res$sets$cs[[cs_iter]]]

				pis = susie_res[['alpha']][cs_num, susie_res$sets$cs[[cs_iter]]]

				posterior_causal_effect_sizes = susie_res[['mu']][cs_num,]*susie_res[['alpha']][cs_num,]

				susie_mu = susie_res[['mu']][cs_num,]

				susie_alpha = susie_res[['alpha']][cs_num,]

				susie_sdev = sqrt(susie_res[['mu2']][cs_num,] - (susie_mu^2))



				lead_variant_index = which(pis==max(pis))[1]
				lead_variant = cs_variant_ids[lead_variant_index]
				lead_variant_position = as.numeric(strsplit(lead_variant, "_")[[1]][2])

				if (lead_variant_position >= window_specific_start_position && lead_variant_position < window_specific_end_position) {
					variant_string = paste0(cs_variant_ids, collapse=",")
					pi_string = paste0(pis, collapse=",")
					new_line <- paste0(chrom_num, "\t", lead_variant_position, "\t", lead_variant, "\t", max(pis), "\t", window_id, "\t", variant_string, "\t", pi_string, "\n")
					cat(new_line)

					# Save posterior mean causal effect sizes for this component
					pmcef_df <- data.frame(variant_id=susie_res$variant_ids, component_posterior_mean_causal_effect_sizes=posterior_causal_effect_sizes, mu=susie_mu, alpha=susie_alpha, mu_sd=susie_sdev)
					pmcef_output_file <- paste0(output_dir, trait_name, "_", lead_variant, "_", window_id, "_component_posterior_mean_causal_effect_sizes.txt")
					write.table(pmcef_df, file=pmcef_output_file, quote=FALSE, sep='\t', row.names = FALSE)
				}
			}
		}
	}
	sink()
}


