args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

tissue_specific_standard_dev_comparison_barplot <- function(ts_var_df, robust_ts_var_df, trait_name) {
	expected_sdev = c(sqrt(1.0/ts_var_df$expected_precision), sqrt(1.0/robust_ts_var_df$expected_precision))
	versions <- c(rep("tgfm", length(ts_var_df$expected_precision)), rep("robust_tgfm", length(robust_ts_var_df$expected_precision)))
	df = data.frame(tissue=c(as.character(ts_var_df$tissue), as.character(robust_ts_var_df$tissue)), std_dev=expected_sdev, version=versions)
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	#df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV−transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c−1="Brain_Spinal_cord")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	df$version = factor(df$version, levels=c("tgfm", "robust_tgfm"))

	p<-ggplot(data=df, aes(x=tissue, y=std_dev, fill=version)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Tissue-specific\nstandard deviation", x="", title=trait_name) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
  	return(p)

}
tissue_specific_standard_dev_comparison_barplot2 <- function(ts_var_df, robust_ts_var_df, trait_name) {
	expected_sdev = c(sqrt(1.0/ts_var_df$expected_precision), sqrt(1.0/robust_ts_var_df$expected_precision))
	versions <- c(rep("robust_tgfm", length(ts_var_df$expected_precision)), rep("robust_mog_tgfm", length(robust_ts_var_df$expected_precision)))
	df = data.frame(tissue=c(as.character(ts_var_df$tissue), as.character(robust_ts_var_df$tissue)), std_dev=expected_sdev, version=versions)
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	#df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV−transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c−1="Brain_Spinal_cord")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	df$version = factor(df$version, levels=c("robust_tgfm", "robust_mog_tgfm"))

	p<-ggplot(data=df, aes(x=tissue, y=std_dev, fill=version)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Tissue-specific\nstandard deviation", x="", title=trait_name) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
  	return(p)

}

tissue_specific_excess_kurtosis_barplot <- function(robust_ts_var_df,trait_name) {
	expected_kurtosis = c(robust_ts_var_df$excess_kurtosis)
	versions <- c(rep("robust_mog_tgfm", length(robust_ts_var_df$excess_kurtosis)))
	df = data.frame(tissue=as.character(robust_ts_var_df$tissue), excess_kurtosis=expected_kurtosis, version=versions)
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	p<-ggplot(data=df, aes(x=tissue, y=excess_kurtosis, fill=version)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Tissue-specific\nkurtosis", x="", title=trait_name) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
  	return(p)
}

make_variance_proportion_heatmap <- function(pseudotissue_gtex_rss_multivariate_twas_dir, trait_names, fusion_weights, robust) {
	variance_proportion_arr <- c()
	tissue_names_arr <- c()
	trait_names_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		if (robust=="robust") {
			precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_", fusion_weights, "_robust_tissue_specific_prior_precision_temp.txt")
		} else {
			precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_", fusion_weights, "_tissue_specific_prior_precision_temp.txt")
		}

		fusion_precision_df <- read.table(precision_file, header=TRUE)

		trait_tissue_names <- as.character(fusion_precision_df$tissue)
		trait_var = 1.0/fusion_precision_df$expected_precision
		norm_trait_var = trait_var/sum(trait_var)
		
		# Add to normalized arr
		variance_proportion_arr <- c(variance_proportion_arr, norm_trait_var)
		tissue_names_arr <- c(tissue_names_arr, trait_tissue_names)
		trait_names_arr <- c(trait_names_arr, rep(trait_name, length(trait_tissue_names)))
	}

	# PUt in inorganized data frame
	df <- data.frame(variance_proportion=variance_proportion_arr, tissue=tissue_names_arr, trait=trait_names_arr)

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	p <- ggplot(df, aes(x = tissue, y = trait, fill = variance_proportion)) +
  		geom_tile() +
  		theme(text = element_text(size=10), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		labs(fill="Proportion of expression-mediated\ntrait heritability",title=paste0("FUSION=", fusion_weights, " / ", robust))
  	return(p)
}

make_comp_nominal_pvalue_thresh_fraction_bar_plot <- function(df, title_string, causal_tissue) {

	df$causal_tissue_bool = as.character(df$tissue) %in% causal_tissue


	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$nominal_fusion_twas_z_score[indices])) > 3.0)

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

	probs_hash <- hash()
	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()
	fraction_std_err_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]

			frac = sum(probs)/length(probs)
			frac_se = std_mean(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)
			fraction_std_err_arr <- c(fraction_std_err_arr, frac_se)

			trait_tissue_name <- paste0(trait_name, "_", tissue_name)
			probs_hash[trait_tissue_name] = probs
		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr, fraction_std_err=fraction_std_err_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
	}

	listy <- list()
	listy[[1]] = df
	listy[[2]] = probs_hash

	return(listy)
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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$tgfm_tissue_prior_twas_z_score[indices])) > 2.0)

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

	probs_hash <- hash()
	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()
	fraction_std_err_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]

			frac = sum(probs)/length(probs)
			frac_se = std_mean(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)
			fraction_std_err_arr <- c(fraction_std_err_arr, frac_se)

			trait_tissue_name <- paste0(trait_name, "_", tissue_name)
			probs_hash[trait_tissue_name] = probs
		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr, fraction_std_err=fraction_std_err_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
	}

	listy <- list()
	listy[[1]] = df
	listy[[2]] = probs_hash

	return(listy)
}
extract_proportion_of_trait_components_mediated_in_each_tissue_robust_tgfm_tissue_prior_multivariate_alpha_z <- function(trait_data_list, trait_names, tissue_names) {
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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$robust_tgfm_tissue_prior_twas_z_score[indices])) > 2.0)

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

	probs_hash <- hash()
	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()
	fraction_std_err_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]

			frac = sum(probs)/length(probs)
			frac_se = std_mean(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)
			fraction_std_err_arr <- c(fraction_std_err_arr, frac_se)

			trait_tissue_name <- paste0(trait_name, "_", tissue_name)
			probs_hash[trait_tissue_name] = probs
		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr, fraction_std_err=fraction_std_err_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
	}

	listy <- list()
	listy[[1]] = df
	listy[[2]] = probs_hash

	return(listy)
}

std_mean <- function(x) sd(x)/sqrt(length(x))

extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior_likelihood_threshold <- function(trait_data_list, trait_names, tissue_names, column_name, likelihood_thresh) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)

	frac_comp <- c()



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]
		
		trait_components = as.character(unique(df$trait_component))
		trait_counter = 0
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(1e-10, num_tiss)
			max_likelihoods <- rep(-1000, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum((component_df[column_name][,1])[indices])
					max_likelihoods[tissue_iter] = max(component_df["robust_tgfm_rss_regression_tissue_prior_log_likelihood"][,1][indices])
				}
			}
			if (max(max_likelihoods) > likelihood_thresh) {
			trait_counter = trait_counter + 1
			tissues_normalized = tissues_un_normalized/sum(tissues_un_normalized)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}
			}
		}
		frac_comp <- c(frac_comp, trait_counter/length(trait_components))
	}

	df_frac <- data.frame(trait=trait_names, fraction_components=frac_comp)
	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	probs_hash <- hash()
	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()
	fraction_std_err_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]

			frac = sum(probs)/length(probs)
			frac_se = std_mean(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)
			fraction_std_err_arr <- c(fraction_std_err_arr, frac_se)

			trait_tissue_name <- paste0(trait_name, "_", tissue_name)
			probs_hash[trait_tissue_name] = probs
		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr, fraction_std_err=fraction_std_err_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
	}

	listy <- list()
	listy[[1]] = df
	listy[[2]] = probs_hash
	listy[[3]] = df_frac$fraction_components


	return(listy)
}

generate_trait_to_loglike_df_hash <- function(trait_data_list, trait_names, tissue_names, column_name) {
	num_tiss <- length(tissue_names)

	log_like_hash <- hash()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df <- trait_data_list[[trait_iter]]
		max_trait_comp_prob <- c()
		max_log_likelihood <- c()
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(1e-10, num_tiss)
			max_likelihoods <- rep(2.3, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum((component_df[column_name][,1])[indices])
					max_likelihoods[tissue_iter] = max(component_df["robust_tgfm_rss_regression_tissue_prior_log_likelihood"][,1][indices])
				}
			
			}

			tissues_normalized = tissues_un_normalized/sum(tissues_un_normalized)

			max_trait_comp_prob <- c(max_trait_comp_prob, max(tissues_normalized))
			max_log_likelihood <- c(max_log_likelihood, max(max_likelihoods))
		}

		df_t <- data.frame(component_prob=max_trait_comp_prob, component_log_like=max_log_likelihood)
		log_like_hash[trait_name] = df_t
	}
	return(log_like_hash)
}

extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior <- function(trait_data_list, trait_names, tissue_names, column_name) {
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
			tissues_un_normalized <- rep(1e-10, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum((component_df[column_name][,1])[indices])

				}
			
			}

			tissues_normalized = tissues_un_normalized/sum(tissues_un_normalized)
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

	probs_hash <- hash()
	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()
	fraction_std_err_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]

			frac = sum(probs)/length(probs)
			frac_se = std_mean(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)
			fraction_std_err_arr <- c(fraction_std_err_arr, frac_se)

			trait_tissue_name <- paste0(trait_name, "_", tissue_name)
			probs_hash[trait_tissue_name] = probs
		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr, fraction_std_err=fraction_std_err_arr)

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
	}

	listy <- list()
	listy[[1]] = df
	listy[[2]] = probs_hash

	return(listy)
}

get_probs_in_causal_tissues_cross_traits <- function(dicti, trait_name, causal_tissue) {
	for (causal_tissue_iter in 1:length(causal_tissue)) {
		causal_tissue_name <- causal_tissue[causal_tissue_iter]
		trait_tissue_name <- paste0(trait_name, "_", causal_tissue_name)
		if (causal_tissue_iter == 1) {
			probs <- dicti[[trait_tissue_name]]
		} else {
			probs <- probs + dicti[[trait_tissue_name]]
		}
	}

	return(probs)
}

make_causal_tissue_prob <- function(methods, prob_list, trait_name, intercept) {
	K <- length(methods)
	method_arr <- c()
	prob_arr <- c()
	se_prob_arr <- c()

	for (kk in 1:K) {
		probs <- prob_list[[kk]]
		method_arr <- c(method_arr, methods[kk])
		prob_arr <- c(prob_arr, mean(probs))
		se_prob_arr <- c(se_prob_arr, std_mean(probs))
	}

	df <- data.frame(method=factor(method_arr, levels=methods), probability=prob_arr, se_probability=se_prob_arr)
	p <- ggplot(df) +
    		geom_bar( aes(x=method, y=probability), stat="identity", fill="skyblue", alpha=0.7) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    		labs(y="Fraction of mediated trait components\nin correct tissue", x="",title=trait_name) +
    		geom_hline(yintercept=intercept, linetype="dashed", color = "black") +
    		geom_errorbar( aes(x=method, ymin=probability-se_probability, ymax=probability+se_probability), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    		figure_theme()
    return(p)
}

make_distribution_of_posterior_probs_in_each_tissue <- function(posterior_probs_hash, trait_name, tissue_names, model_name, causal_tissue) {
	tissue_names <- as.character(tissue_names)

	tissue_arr <- c()
	posterior_prob_arr <- c()
	causal_tissue_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]

		probs <- posterior_probs_hash[[paste0(trait_name, "_", tissue_name)]]

		posterior_prob_arr <- c(posterior_prob_arr, probs)
		tissue_arr <- c(tissue_arr, rep(tissue_name, length(probs)))
	}

	df <- data.frame(tissue=tissue_arr, posterior_prob=posterior_prob_arr)
	p <- ggplot(df, aes(x=posterior_prob)) + geom_histogram() + geom_histogram() + figure_theme() +
		facet_wrap( ~ tissue, ncol = 4)
	#p <- ggplot(df, aes(x=tissue, y=posterior_prob, fill=causal_tissue)) +
	#	figure_theme() + 
	#	geom_boxplot() +
	#	theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
	return(p)
}

create_component_log_like_histogram <- function(df, trait_name) {
	p <- ggplot(df, aes(x=component_log_like)) + geom_histogram() + figure_theme() +
	labs(title=trait_name, x="Component max log likelihood")

	return(p)
}

create_component_log_like_causal_prob_scatter <- function(df, trait_name) {
	df$component_log_like[df$component_log_like > 100] = 100.0
	p <- ggplot(df, aes(x=component_log_like, y=component_prob)) + geom_point() + figure_theme() +
	labs(title=trait_name, x="Component max log likelihood", y="Component max probability")

	return(p)
}


extract_causal_tissue_prob <- function(tissue_df, trait_name, causal_tissue) {
	causal_prob = 0.0
	trait_tissue_df = tissue_df[as.character(tissue_df$trait)==trait_name,]
	for (tissue_iter in 1:length(trait_tissue_df$tissue)) {
		if (trait_tissue_df$tissue[tissue_iter] %in% causal_tissue) {
			causal_prob = causal_prob + trait_tissue_df$fraction_nominal_twas_associations[tissue_iter]
		}
	}
	return(causal_prob)

}

fraction_components_passed_stacked_barplot <- function(fraction_components_na, fraction_components_1, fraction_components_3, fraction_components_10, fraction_components_20, trait_name) {
	fractions <- c()
	methods <- c()
	correct_tissue <- c()

	fractions <- c(fractions, fraction_components_na, 1.0-fraction_components_na)
	fractions <- c(fractions, fraction_components_1, 1.0-fraction_components_1)
	fractions <- c(fractions, fraction_components_3, 1.0-fraction_components_3)
	fractions <- c(fractions, fraction_components_10,  1.0-fraction_components_10)
	fractions <- c(fractions, fraction_components_20, 1.0-fraction_components_20)
	methods <- c(methods, "log_like_na", "log_like_na")
	methods <- c(methods, "log_like_1", "log_like_1")
	methods <- c(methods, "log_like_3", "log_like_3")
	methods <- c(methods, "log_like_10", "log_like_10")
	methods <- c(methods, "log_like_20", "log_like_20")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")


	df <- data.frame(fraction=fractions, method=factor(methods, levels=c("log_like_na","log_like_1", "log_like_3","log_like_10","log_like_20")), correct_tissue=correct_tissue)


	p <- ggplot(df, aes(x = method, y = fraction, fill = correct_tissue)) + 
 		 geom_bar(stat = "identity") +
 		 scale_fill_manual(values = c("#DADAEB", "#9E9AC8")) +
  		figure_theme() +
  		labs(y="Fraction components", x="",title=trait_name, fill="passed_filter")


  	return(p)	
}

fraction_correct_stacked_barplot <- function(causal_tissue_prob_na, causal_tissue_prob_1, causal_tissue_prob_3, causal_tissue_prob_10, causal_tissue_prob_20, trait_name) {
	fractions <- c()
	methods <- c()
	correct_tissue <- c()

	fractions <- c(fractions, causal_tissue_prob_na, 1.0-causal_tissue_prob_na)
	fractions <- c(fractions, causal_tissue_prob_1, 1.0-causal_tissue_prob_1)
	fractions <- c(fractions, causal_tissue_prob_3, 1.0-causal_tissue_prob_3)
	fractions <- c(fractions, causal_tissue_prob_10, 1.0-causal_tissue_prob_10)
	fractions <- c(fractions, causal_tissue_prob_20, 1.0-causal_tissue_prob_20)
	methods <- c(methods, "log_like_na", "log_like_na")
	methods <- c(methods, "log_like_1", "log_like_1")
	methods <- c(methods, "log_like_3", "log_like_3")
	methods <- c(methods, "log_like_10", "log_like_10")
	methods <- c(methods, "log_like_20", "log_like_20")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")
	correct_tissue <- c(correct_tissue, "True", "False")

	df <- data.frame(fraction=fractions, method=factor(methods, levels=c("log_like_na","log_like_1","log_like_3","log_like_10","log_like_20")), correct_tissue=correct_tissue)


	p <- ggplot(df, aes(x = method, y = fraction, fill = correct_tissue)) + 
 		 geom_bar(stat = "identity") +
  		scale_fill_manual(values = c("#DADAEB", "#9E9AC8")) +
  		figure_theme() +
  		labs(y="Fraction components", x="",title=trait_name)


  	return(p)
}


trait_file <- args[1]
gtex_pseudotissue_file <- args[2]
pseudotissue_gtex_rss_multivariate_twas_dir <- args[3]
rss_multivariate_twas_visualization_dir <- args[4]

# Load in gtex tissues
tissue_df <- read.table(gtex_pseudotissue_file,header=TRUE)
tissue_names <- as.character(tissue_df$pseudotissue_name)
tissue_names <- str_replace_all(tissue_names, "-", "_")

# trait names
trait_df <- read.table(trait_file, header=TRUE)

trait_df = trait_df[as.character(trait_df$study_name) != "biochemistry_Triglycerides", ]
trait_names <- as.character(trait_df$study_name)
trait_names <- c(trait_names[3:4], trait_names[8], trait_names[10])

# Visualize variance componenets with bar plot
if (FALSE) {
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	# Data files
	fusion_true_precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_True_tissue_specific_prior_precision_temp.txt")
	fusion_true_robust_precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_True_robust_tissue_specific_prior_precision_temp.txt")
	fusion_false_precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_False_tissue_specific_prior_precision_temp.txt")
	fusion_false_robust_precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_False_robust_tissue_specific_prior_precision_temp.txt")

	fusion_true_precision_df <- read.table(fusion_true_precision_file, header=TRUE)
	fusion_true_robust_precision_df <- read.table(fusion_true_robust_precision_file, header=TRUE)
	fusion_false_precision_df <- read.table(fusion_false_precision_file, header=TRUE)
	fusion_false_robust_precision_df <- read.table(fusion_false_robust_precision_file, header=TRUE)

	bar_plot_fusion_false = tissue_specific_standard_dev_comparison_barplot(fusion_false_precision_df, fusion_false_robust_precision_df, paste0(trait_name, " / SuSiE distribution eQTL"))
	bar_plot_fusion_true = tissue_specific_standard_dev_comparison_barplot(fusion_true_precision_df, fusion_true_robust_precision_df, paste0(trait_name, " / FUSION eQTL point estimates"))

	bar_plot <- plot_grid(bar_plot_fusion_true, bar_plot_fusion_false, nrow=2)

	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name , "_standard_vs_robust_tissue_specific_variance_barplot.pdf")
	ggsave(bar_plot, file=output_file, width=7.2, height=7.5, units="in")
}
}


# Visualize proportion of  variance componenets coming from each tissue with heatmap
if (FALSE) {
fusion_weights="False"
robust="robust"
output_file <- paste0(rss_multivariate_twas_visualization_dir , "tissue_variance_proportion_fusion_", fusion_weights, "_", robust, ".pdf")
heatmap <- make_variance_proportion_heatmap(pseudotissue_gtex_rss_multivariate_twas_dir, trait_names, fusion_weights, robust)
ggsave(heatmap, file=output_file, width=7.2, height=5.0, units="in")

fusion_weights="False"
robust="non_robust"
output_file <- paste0(rss_multivariate_twas_visualization_dir , "tissue_variance_proportion_fusion_", fusion_weights, "_", robust, ".pdf")
heatmap <- make_variance_proportion_heatmap(pseudotissue_gtex_rss_multivariate_twas_dir, trait_names, fusion_weights, robust)
ggsave(heatmap, file=output_file, width=7.2, height=5.0, units="in")

fusion_weights="True"
robust="robust"
output_file <- paste0(rss_multivariate_twas_visualization_dir , "tissue_variance_proportion_fusion_", fusion_weights, "_", robust, ".pdf")
heatmap <- make_variance_proportion_heatmap(pseudotissue_gtex_rss_multivariate_twas_dir, trait_names, fusion_weights, robust)
ggsave(heatmap, file=output_file, width=7.2, height=5.0, units="in")

fusion_weights="True"
robust="non_robust"
output_file <- paste0(rss_multivariate_twas_visualization_dir , "tissue_variance_proportion_fusion_", fusion_weights, "_", robust, ".pdf")
heatmap <- make_variance_proportion_heatmap(pseudotissue_gtex_rss_multivariate_twas_dir, trait_names, fusion_weights, robust)
ggsave(heatmap, file=output_file, width=7.2, height=5.0, units="in")
}


#trait_names <- c("bp_DIASTOLICadjMEDz")

fusion_weights="False"

trait_to_tgfm_ch_data_hash = list()
for (trait_iter in 1:length(trait_names)) {
	gene_version = "cis_heritable_genes"
	trait_name <- trait_names[trait_iter]
	trait_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_component_organized_tgfm_results.txt")
	df <- read.table(trait_file, header=TRUE)
	trait_to_tgfm_ch_data_hash[[trait_iter]] = df
}



# Distribution of posterior probabilities per tissue histogram
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob")
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_probs_hash = proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list[[2]]

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	p_robust_tgfm <- make_distribution_of_posterior_probs_in_each_tissue(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_probs_hash, trait_name, tissue_names, "Robust_TGFM")
	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_robust_tgfm_posterior_prob_distribution_per_tissue_histogram.pdf")
	ggsave(p_robust_tgfm, file=output_file, width=7.2, height=8.5, units="in")
}



trait_to_loglike_df_hash <- generate_trait_to_loglike_df_hash(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob")


for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	log_like_df <- trait_to_loglike_df_hash[[trait_name]]

	component_log_like_histo = create_component_log_like_histogram(log_like_df, trait_name)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_max_component_log_likelihood_histogram.pdf")
	ggsave(component_log_like_histo, file=output_file, width=7.2, height=3.5, units="in")

	component_log_like_histo = create_component_log_like_causal_prob_scatter(log_like_df, trait_name)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_max_component_log_likelihood_causal_prob.pdf")
	ggsave(component_log_like_histo, file=output_file, width=7.2, height=5.5, units="in")

}


proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_1 = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior_likelihood_threshold(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob", 1.0)
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_3 = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior_likelihood_threshold(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob", 3.0)
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_10 = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior_likelihood_threshold(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob", 10.0)
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_20 = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior_likelihood_threshold(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob", 20.0)



for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	if (trait_name == "blood_WHITE_COUNT") {
		causal_tissue =c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "body_WHRadjBMIz") {
		causal_tissue = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum")
	}
	if (trait_name == "biochemistry_Cholesterol") {
		causal_tissue = c("Liver")
	}
	if (trait_name == "bp_DIASTOLICadjMEDz") {
		causal_tissue = c("Artery_Heart", "Artery_Tibial", "Heart_Atrial_Appendage", "Heart_Left_Ventricle")
	}
	if (trait_name == "blood_EOSINOPHIL_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "blood_RBC_DISTRIB_WIDTH") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "blood_RED_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "cov_EDU_COLLEGE") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
	}
	if (trait_name == "mental_NEUROTICISM") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
	}
	if (trait_name == "Disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}

	causal_tissue_prob_na = extract_causal_tissue_prob(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list[[1]], trait_name, causal_tissue)
	causal_tissue_prob_1 = extract_causal_tissue_prob(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_1[[1]], trait_name, causal_tissue)
	causal_tissue_prob_3 = extract_causal_tissue_prob(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_3[[1]], trait_name, causal_tissue)
	causal_tissue_prob_10 = extract_causal_tissue_prob(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_10[[1]], trait_name, causal_tissue)
	causal_tissue_prob_20 = extract_causal_tissue_prob(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_20[[1]], trait_name, causal_tissue)

	fraction_components_na = 1.0
	fraction_components_1 = proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_1[[3]][trait_iter]
	fraction_components_3 = proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_3[[3]][trait_iter]
	fraction_components_10 = proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_10[[3]][trait_iter]
	fraction_components_20 = proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list_thresh_20[[3]][trait_iter]

	frac_correct_barplot <- fraction_correct_stacked_barplot(causal_tissue_prob_na, causal_tissue_prob_1, causal_tissue_prob_3, causal_tissue_prob_10, causal_tissue_prob_20, trait_name)
	frac_component_barplot <- fraction_components_passed_stacked_barplot(fraction_components_na, fraction_components_1, fraction_components_3,fraction_components_10, fraction_components_20, trait_name)

	p_merged = plot_grid(frac_component_barplot,frac_correct_barplot,ncol=1)


	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_robust_tgfm_fraction_correct_barplot_thresholds.pdf")
	ggsave(p_merged, file=output_file, width=7.2, height=7.5, units="in")


}





if (FALSE) {
proportion_trait_components_mediated_tgfm_univariate_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_univariate_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('one')
proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_multivariate_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('two')
proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_robust_tgfm_tissue_prior_multivariate_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('three')
proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "tgfm_rss_regression_tissue_prior_posterior_prob")
print('four')
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob")
print('five')




proportion_trait_components_mediated_tgfm_univariate_alpha_z = proportion_trait_components_mediated_tgfm_univariate_alpha_z_list[[1]]
proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z = proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z_list[[1]]
proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z = proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z_list[[1]]
proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior = proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior_list[[1]]
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior = proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list[[1]]




for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	if (trait_name == "blood_WHITE_COUNT") {
		causal_tissue =c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "body_WHRadjBMIz") {
		causal_tissue = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum")
	}
	if (trait_name == "biochemistry_Cholesterol") {
		causal_tissue = c("Liver")
	}
	if (trait_name == "bp_DIASTOLICadjMEDz") {
		causal_tissue = c("Artery_Heart", "Artery_Tibial", "Heart_Atrial_Appendage", "Heart_Left_Ventricle")
	}
	if (trait_name == "blood_EOSINOPHIL_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "blood_RBC_DISTRIB_WIDTH") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "blood_RED_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "cov_EDU_COLLEGE") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
	}
	if (trait_name == "mental_NEUROTICISM") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
	}
	if (trait_name == "Disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}

	p_univariate <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_univariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_univariate_alpha_z$trait)==trait_name,], "univariate_twas", causal_tissue)
	p_multivariate_tp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z$trait)==trait_name,], "multivariate_tissue_specific_prior_twas", causal_tissue)
	p_robust_multivariate_tp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z$trait)==trait_name,], "multivariate_tissue_specific_prior_robust_twas", causal_tissue)

	p_tgfm <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior[as.character(proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior$trait)==trait_name,], "TGFM", causal_tissue)
	p_robust_tgfm <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior[as.character(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior$trait)==trait_name,], "Robust_TGFM", causal_tissue)

	p_merged = plot_grid(p_univariate,NULL, p_multivariate_tp, p_robust_multivariate_tp, p_tgfm,p_robust_tgfm,ncol=2)

	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_fraction_causal_mediated.pdf")
	ggsave(p_merged, file=output_file, width=10.7, height=9.5, units="in")
}




for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	methods <- c("univariate_twas", "multivariate_tissue_specific_prior_twas", "multivariate_tissue_specific_prior_robust_twas", "TGFM", "Robust_TGFM")

	if (trait_name == "blood_WHITE_COUNT") {
		causal_tissue =c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "body_WHRadjBMIz") {
		causal_tissue = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum")
	}
	if (trait_name == "biochemistry_Cholesterol") {
		causal_tissue = c("Liver")
	}
	if (trait_name == "bp_DIASTOLICadjMEDz") {
		causal_tissue = c("Artery_Heart", "Artery_Tibial", "Heart_Atrial_Appendage", "Heart_Left_Ventricle")
	}
	if (trait_name == "blood_EOSINOPHIL_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "blood_RBC_DISTRIB_WIDTH") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "blood_RED_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}
	if (trait_name == "cov_EDU_COLLEGE") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
	}
	if (trait_name == "mental_NEUROTICISM") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
	}
	if (trait_name == "Disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}


	prob_list <- list()
	prob_list[[1]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_univariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[2]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[3]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[4]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior_list[[2]], trait_name, causal_tissue)
	prob_list[[5]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list[[2]], trait_name, causal_tissue)


	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_fraction_causal_mediated_in_correct_tiss.pdf")
	causal_tissue_prob_se_plot <- make_causal_tissue_prob(methods, prob_list, trait_name, length(causal_tissue)/38.0)
	ggsave(causal_tissue_prob_se_plot, file=output_file, width=7.2, height=5.0, units="in")

}





}

























































































############################
# Excess kurtosis stuff
##########################
if (FALSE) {
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	# Data files
	fusion_false_mog_robust_precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_False_robust_mog_tissue_specific_prior_precision_temp.txt")
	fusion_false_robust_precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_False_robust_tissue_specific_prior_precision_temp.txt")


	fusion_false_robust_mog_precision_df <- read.table(fusion_false_mog_robust_precision_file, header=TRUE)
	fusion_false_robust_precision_df <- read.table(fusion_false_robust_precision_file, header=TRUE)

	bar_plot_fusion_false = tissue_specific_standard_dev_comparison_barplot2(fusion_false_robust_precision_df, fusion_false_robust_mog_precision_df, paste0(trait_name, " / SuSiE distribution eQTL"))

	bar_plot_excess_kurtosis = tissue_specific_excess_kurtosis_barplot(fusion_false_robust_mog_precision_df, paste0(trait_name))

	bar_plot <- plot_grid(bar_plot_fusion_false, bar_plot_excess_kurtosis, nrow=2)

	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name , "_robust_vs_robust_mog_tissue_specific_variance_barplot.pdf")
	ggsave(bar_plot, file=output_file, width=7.2, height=7.5, units="in")
}
}

