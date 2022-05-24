args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(reshape2)
library(pheatmap)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_univariate_alpha_z <- function(trait_data_list, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]

		tissue_counts <- rep(0, num_tiss)
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- tissue_names[tissue_iter]
			num_genes <- sum(df$tissue==tissue_name)
			tissue_counts[tissue_iter] = num_genes
		}
		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$nominal_tgfm_twas_z_score[indices])) > 3.0)

				}
			
			}

			tissues_un_normalized = tissues_un_normalized + pseudocount
			tissues_normalized =tissues_un_normalized/sum(tissues_un_normalized)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(df[as.character(df$trait)==trait_name,])
	}

	return(df)
}

extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_constant_prior_multivariate_alpha_z <- function(trait_data_list, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]

		tissue_counts <- rep(0, num_tiss)
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- tissue_names[tissue_iter]
			num_genes <- sum(df$tissue==tissue_name)
			tissue_counts[tissue_iter] = num_genes
		}
		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$multivariate_tgfm_const_prior_twas_z_score[indices])) > 2.0)

				}
			
			}

			tissues_un_normalized = tissues_un_normalized + pseudocount
			tissues_normalized =tissues_un_normalized/sum(tissues_un_normalized)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(df[as.character(df$trait)==trait_name,])
	}

	return(df)
}


extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_multivariate_alpha_z <- function(trait_data_list, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]

		tissue_counts <- rep(0, num_tiss)
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- tissue_names[tissue_iter]
			num_genes <- sum(df$tissue==tissue_name)
			tissue_counts[tissue_iter] = num_genes
		}
		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$multivariate_tgfm_tissue_prior_twas_z_score[indices])) > 2.0)

				}
			
			}

			tissues_un_normalized = tissues_un_normalized + pseudocount
			tissues_normalized =tissues_un_normalized/sum(tissues_un_normalized)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(df[as.character(df$trait)==trait_name,])
	}

	return(df)
}




extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_sum_posterior_thresholded <- function(trait_data_list, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]
		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(component_df$tgfm_const_prior_max_posterior_prob[indices]) > .5)

				}
			
			}

			tissues_un_normalized = tissues_un_normalized + pseudocount
			tissues_normalized =tissues_un_normalized/sum(tissues_un_normalized)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(df[as.character(df$trait)==trait_name,])
	}

	return(df)
}


extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_sum_posterior <- function(trait_data_list, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]
		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum(component_df$tgfm_tissue_prior_sum_posterior_prob[indices])

				}
			
			}

			tissues_normalized = tissues_un_normalized
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(df[as.character(df$trait)==trait_name,])
	}

	return(df)
}

extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_constant_prior_sum_posterior <- function(trait_data_list, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]
		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum(component_df$tgfm_const_prior_sum_posterior_prob[indices])

				}
			
			}

			tissues_normalized = tissues_un_normalized
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(df[as.character(df$trait)==trait_name,])
	}

	return(df)
}


nominal_fusion_z_score_vs_nominal_tgfm_z_score <- function(nominal_fusion_twas_z_score, nominal_tgfm_twas_z_score, trait_name) {
	nominal_fusion_twas_z_score[nominal_fusion_twas_z_score < -40] = -40
	nominal_fusion_twas_z_score[nominal_fusion_twas_z_score > 40] = 40
	nominal_tgfm_twas_z_score[nominal_tgfm_twas_z_score < -40] = -40
	nominal_tgfm_twas_z_score[nominal_tgfm_twas_z_score > 40] = 40

	df <- data.frame(nominal_fusion_twas_z_score=nominal_fusion_twas_z_score, nominal_tgfm_twas_z_score=nominal_tgfm_twas_z_score)

	corry <- cor(df$nominal_fusion_twas_z_score, df$nominal_tgfm_twas_z_score)
	p <- ggplot(df, aes(x=nominal_fusion_twas_z_score, y=nominal_tgfm_twas_z_score)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline() +
		labs(title=paste0(trait_name, " / pearson correlation: ", corry))
	return(p)	
}

nominal_fusion_z_score_vs_nominal_tgfm_z_score_colored_by_stochastic_sd <- function(nominal_fusion_twas_z_score, nominal_tgfm_twas_z_score, stochastic_twas_z_score_sd, trait_name) {
	nominal_fusion_twas_z_score[nominal_fusion_twas_z_score < -40] = -40
	nominal_fusion_twas_z_score[nominal_fusion_twas_z_score > 40] = 40
	nominal_tgfm_twas_z_score[nominal_tgfm_twas_z_score < -40] = -40
	nominal_tgfm_twas_z_score[nominal_tgfm_twas_z_score > 40] = 40

	df <- data.frame(nominal_fusion_twas_z_score=nominal_fusion_twas_z_score, nominal_tgfm_twas_z_score=nominal_tgfm_twas_z_score, stochastic_twas_z_score_sd=stochastic_twas_z_score_sd)
	df = df[!is.na(df$stochastic_twas_z_score_sd),]
	
	df$stochastic_twas_z_score_sd[df$stochastic_twas_z_score_sd > 40] =40

	#df$stochastic_twas_z_score_sd = log(abs(df$nominal_fusion_twas_z_score/df$stochastic_twas_z_score_sd))

	corry <- cor(df$nominal_fusion_twas_z_score, df$nominal_tgfm_twas_z_score)
	p <- ggplot(df, aes(x=nominal_fusion_twas_z_score, y=nominal_tgfm_twas_z_score, color=stochastic_twas_z_score_sd)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline() +
		labs(title=paste0(trait_name, " / pearson correlation: ", corry))
	return(p)		
}
nominal_fusion_z_score_vs_nominal_tgfm_z_score_colored_by_gene_sd <- function(nominal_fusion_twas_z_score, nominal_tgfm_twas_z_score, gene_sd, trait_name) {
	nominal_fusion_twas_z_score[nominal_fusion_twas_z_score < -40] = -40
	nominal_fusion_twas_z_score[nominal_fusion_twas_z_score > 40] = 40
	nominal_tgfm_twas_z_score[nominal_tgfm_twas_z_score < -40] = -40
	nominal_tgfm_twas_z_score[nominal_tgfm_twas_z_score > 40] = 40

	df <- data.frame(nominal_fusion_twas_z_score=nominal_fusion_twas_z_score, nominal_tgfm_twas_z_score=nominal_tgfm_twas_z_score, abs_diff=abs(nominal_fusion_twas_z_score-nominal_tgfm_twas_z_score), gene_sd=gene_sd)
	

	#df$stochastic_twas_z_score_sd = log(abs(df$nominal_fusion_twas_z_score/df$stochastic_twas_z_score_sd))

	corry <- cor(df$nominal_fusion_twas_z_score, df$nominal_tgfm_twas_z_score)
	#p <- ggplot(df, aes(x=nominal_fusion_twas_z_score, y=nominal_tgfm_twas_z_score, color=gene_sd)) + geom_point(size=.1) +
	p <- ggplot(df, aes(x=abs_diff, y=gene_sd, color=gene_sd)) + geom_point(size=.1) +
		figure_theme() +
		#geom_abline() +
		labs(title=paste0(trait_name, " / pearson correlation: ", corry))
	return(p)		
}

multivaraite_tgfm_z_score_vs_nominal_tgfm_z_score <- function(multivariate_tgfm_const_prior_twas_z_score, nominal_tgfm_twas_z_score, trait_name) {
	multivariate_tgfm_const_prior_twas_z_score[multivariate_tgfm_const_prior_twas_z_score > 20] = 20
	multivariate_tgfm_const_prior_twas_z_score[multivariate_tgfm_const_prior_twas_z_score < -20] = -20
	df <- data.frame(multivariate_tgfm_twas_z_score=multivariate_tgfm_const_prior_twas_z_score, nominal_tgfm_twas_z_score=nominal_tgfm_twas_z_score)
	corry <- cor(df$multivariate_tgfm_twas_z_score, df$nominal_tgfm_twas_z_score)
	p <- ggplot(df, aes(x=nominal_tgfm_twas_z_score, y=multivariate_tgfm_twas_z_score)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline(size=.1) +
		geom_hline(yintercept=0,size=.1) + 
		geom_vline(xintercept=0,size=.1) +
		labs(title=paste0(trait_name, " / pearson correlation: ", corry))
	return(p)
}


multivaraite_tgfm_z_score_vs_multivariate_tgfm_z_score_tissue_prior <- function(multivariate_tgfm_const_prior_twas_z_score, multivariate_tgfm_tissue_prior_twas_z_score, tissue_indicator, trait_name, tissue_name) {
	multivariate_tgfm_const_prior_twas_z_score[multivariate_tgfm_const_prior_twas_z_score > 30] = 30
	multivariate_tgfm_const_prior_twas_z_score[multivariate_tgfm_const_prior_twas_z_score < -30] = -30
	multivariate_tgfm_tissue_prior_twas_z_score[multivariate_tgfm_tissue_prior_twas_z_score > 30] = 30
	multivariate_tgfm_tissue_prior_twas_z_score[multivariate_tgfm_tissue_prior_twas_z_score < -30] = -30

	df <- data.frame(multivariate_tgfm_const_prior_twas_z_score=multivariate_tgfm_const_prior_twas_z_score, multivariate_tgfm_tissue_prior_twas_z_score=multivariate_tgfm_tissue_prior_twas_z_score, tissue_indicator=factor(tissue_indicator,levels=c(TRUE,FALSE)))
	corry <- cor(df$multivariate_tgfm_const_prior_twas_z_score, df$multivariate_tgfm_tissue_prior_twas_z_score)
	p <- ggplot(df, aes(x=multivariate_tgfm_const_prior_twas_z_score, y=multivariate_tgfm_tissue_prior_twas_z_score, color=tissue_indicator)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline(size=.1) +
		geom_hline(yintercept=0,size=.1) + 
		geom_vline(xintercept=0,size=.1) +
		labs(title=paste0(trait_name, " / pearson correlation: ", corry), color=paste0(tissue_name, " indicator"), x="Multivariate TGFM TWAS Z score (constant prior)", y="Multivariate TGFM TWAS Z score (tissue prior)")
	return(p)
}


tissue_specific_variance_barplot <- function(ts_var_df, trait_name) {
	expected_var = 1.0/ts_var_df$expected_precision
	df = data.frame(tissue=ts_var_df$tissue, variance=expected_var)
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn")

	p<-ggplot(data=df, aes(x=tissue, y=variance)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Tissue-specific variance", x="", title=trait_name) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
  	return(p)
}
make_comp_nominal_pvalue_thresh_fraction_bar_plot <- function(df, title_string, causal_tissue) {

	df$causal_tissue_bool = as.character(df$tissue) == causal_tissue

	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn")

	indices <- order(df$fraction_nominal_twas_associations)

	tissue_order = as.character(df$tissue)[indices]

	df$tissue = factor(df$tissue, levels=tissue_order)



	p<-ggplot(data=df, aes(x=tissue, y=fraction_nominal_twas_associations, fill=causal_tissue_bool)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="none") +
  			labs(y="Fraction of mediated trait components", x="",title=title_string) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))

  	return(p)

}

gtex_tissue_file <- args[1]
gtex_rss_twas_associations_dir <- args[2]
visualize_twas_results_dir <- args[3]



# Load in gtex tissues
tissue_df <- read.table(gtex_tissue_file,header=TRUE)
tissue_names <- as.character(tissue_df$pseudotissue_name)

# names of traits 
trait_names <- c("blood_WHITE_COUNT", "lung_FEV1FVCzSMOKE", "body_WHRadjBMIz", "bp_DIASTOLICadjMEDz")



for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	tissue_specific_variance_file = paste0(gtex_rss_twas_associations_dir, trait_name, "_tissue_specific_prior_precision_v2_temp.txt")
	ts_var_df = read.table(tissue_specific_variance_file, header=TRUE)
	bar_plot = tissue_specific_variance_barplot(ts_var_df, trait_name)
	output_file <- paste0(visualize_twas_results_dir, trait_name,"_tissue_specific_variance_barplot.pdf")
	ggsave(bar_plot, file=output_file, width=7.2, height=3.5, units="in")


}



# Create dictionary mapping from trait name to tgfm data frame
trait_to_tgfm_data_hash = list()
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_file <- paste0(gtex_rss_twas_associations_dir, trait_name, "_component_organized_tgfm_results.txt")
	df <- read.table(trait_file, header=TRUE)
	trait_to_tgfm_data_hash[[trait_iter]] = df
}


for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_df <- trait_to_tgfm_data_hash[[trait_iter]]

	scatter <- nominal_fusion_z_score_vs_nominal_tgfm_z_score(trait_df$nominal_fusion_twas_z_score, trait_df$nominal_tgfm_twas_z_score, trait_name)
	output_file <- paste0(visualize_twas_results_dir, trait_name,"_nominal_fusion_twas_z_score_vs_nominal_tgfm_twas_z_score_scatter.pdf")
	ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")


	scatter <- multivaraite_tgfm_z_score_vs_nominal_tgfm_z_score(trait_df$multivariate_tgfm_const_prior_twas_z_score, trait_df$nominal_tgfm_twas_z_score, trait_name)
	output_file <- paste0(visualize_twas_results_dir, trait_name,"_multivariate_tgfm_twas_z_score_vs_nominal_tgfm_twas_z_score_scatter.pdf")
	ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")

	if (trait_name == "body_WHRadjBMIz") {
		tissue_name="Adipose_Subcutaneous"
		scatter <- multivaraite_tgfm_z_score_vs_multivariate_tgfm_z_score_tissue_prior(trait_df$multivariate_tgfm_const_prior_twas_z_score, trait_df$multivariate_tgfm_tissue_prior_twas_z_score,trait_df$tissue==tissue_name,trait_name, tissue_name)
		output_file <- paste0(visualize_twas_results_dir, trait_name,"_multivariate_tgfm_twas_z_score_const_prior_vs_multivariate_tgfm_twas_z_score_tissue_prior_scatter_colored_by_", tissue_name, ".pdf")
		ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")
	}


}


proportion_trait_components_mediated_tgfm_univariate_alpha_z = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_univariate_alpha_z(trait_to_tgfm_data_hash, trait_names, tissue_names)
proportion_trait_components_mediated_tgfm_const_prior_multivariate_alpha_z = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_constant_prior_multivariate_alpha_z(trait_to_tgfm_data_hash, trait_names, tissue_names)
proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_multivariate_alpha_z(trait_to_tgfm_data_hash, trait_names, tissue_names)
proportion_trait_components_mediated_tgfm_tissue_prior_sum_posterior = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_sum_posterior(trait_to_tgfm_data_hash, trait_names, tissue_names)


for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	if (trait_name == "blood_WHITE_COUNT") {
		causal_tissue ="Whole_Blood"
	}
	if (trait_name == "body_WHRadjBMIz") {
		causal_tissue = "Adipose_Subcutaneous"
	}
	if (trait_name == "lung_FEV1FVCzSMOKE") {
		causal_tissue = "Esophagus_Muscularis"
	}
	if (trait_name == "bp_DIASTOLICadjMEDz") {
		causal_tissue = "Artery_Tibial"
	}

	p_univariate <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_univariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_univariate_alpha_z$trait)==trait_name,], "univariate_twas", causal_tissue)
	p_multivariate <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_const_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_const_prior_multivariate_alpha_z$trait)==trait_name,], "multivariate_twas", causal_tissue)
	p_multivariate_tp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z$trait)==trait_name,], "multivariate_tissue_specific_prior_twas", causal_tissue)
	p_tgfm <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_tissue_prior_sum_posterior[as.character(proportion_trait_components_mediated_tgfm_tissue_prior_sum_posterior$trait)==trait_name,], "TGFM", causal_tissue)

	p_merged = plot_grid(p_univariate, p_multivariate, p_multivariate_tp, p_tgfm,ncol=2)

	output_file <- paste0(visualize_twas_results_dir, trait_name,"_fraction_causal_mediated.pdf")
	ggsave(p_merged, file=output_file, width=10.7, height=7.5, units="in")

}



#proportion_trait_components_mediated_tgfm_constant_prior_sum_posterior = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_constant_prior_sum_posterior(trait_to_tgfm_data_hash, trait_names, tissue_names)
#proportion_trait_components_mediated_tgfm_tissue_prior_sum_posterior = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_sum_posterior(trait_to_tgfm_data_hash, trait_names, tissue_names)
#proportion_trait_components_mediated_tgfm_tissue_prior_sum_posterior_thresh = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_sum_posterior_thresholded(trait_to_tgfm_data_hash, trait_names, tissue_names)





# Competitive DFs for:
# nominal_fusion_twas_z_score
# nominal_tgfm_twas_z_score
# multivariate_tgfm_const_prior_twas_z_score
# multivariate_tgfm_tissue_prior_twas_z_score
# tgfm_const_prior_sum_posterior_prob
# tgfm_const_prior_max_posterior_prob
# tgfm_tissue_prior_sum_posterior_probtgfm_tissue_prior_max_posterior_prob