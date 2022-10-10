args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
options(warn=1)

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


make_average_causal_prob_for_variant_components_with_nearby_gene_components_bar_plot <- function(df_causal_prob, title_string, causal_tissue, tissue_names) {
	df_causal_prob <- df_causal_prob[df_causal_prob$gene_component_boolean==1.0,]
	tissue_arr <- c()
	average_causal_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		avg_prob <- sum(df_causal_prob$tissue_causal_prob[as.character(df_causal_prob$tissue)==tissue_name])/sum(as.character(df_causal_prob$tissue)==tissue_name)
		tissue_arr <- c(tissue_arr, tissue_name)
		average_causal_arr <- c(average_causal_arr, avg_prob)
	}

	df <- data.frame(tissue=tissue_arr, average_causal_prob=average_causal_arr)

	df$causal_tissue_bool = as.character(df$tissue) %in% causal_tissue

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	indices <- order(df$average_causal_prob)

	tissue_order = as.character(df$tissue)[indices]

	df$tissue = factor(df$tissue, levels=tissue_order)
	
	p<-ggplot(data=df, aes(x=tissue, y=average_causal_prob, fill=causal_tissue_bool)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="none") +
  			labs(y="Fraction of\n trait components\nw gene components", x="",title=title_string) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
  	return(p)
}

make_average_causal_prob_bar_plot <- function(df_causal_prob, title_string, causal_tissue, tissue_names) {
	tissue_arr <- c()
	average_causal_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		avg_prob <- sum(df_causal_prob$tissue_causal_prob[as.character(df_causal_prob$tissue)==tissue_name])/sum(as.character(df_causal_prob$tissue)==tissue_name)
		tissue_arr <- c(tissue_arr, tissue_name)
		average_causal_arr <- c(average_causal_arr, avg_prob)
	}

	df <- data.frame(tissue=tissue_arr, average_causal_prob=average_causal_arr)

	df$causal_tissue_bool = as.character(df$tissue) %in% causal_tissue

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	indices <- order(df$average_causal_prob)

	tissue_order = as.character(df$tissue)[indices]

	df$tissue = factor(df$tissue, levels=tissue_order)
	
	p<-ggplot(data=df, aes(x=tissue, y=average_causal_prob, fill=causal_tissue_bool)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="none") +
  			labs(y="Fraction of\n trait components", x="",title=title_string) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
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
  			labs(y="Fraction of\n trait components", x="",title=title_string) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))

  	return(p)

}

extract_univariate_z_df <- function(trait_data_list, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	raw_arr <- c()
	observed_arr <- c()


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
			tissue_z <- rep(0.0, num_tiss)
			tissue_observed <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$nominal_fusion_twas_z_score[indices])) > 3.0)
					tissue_z[tissue_iter] = max(abs(component_df$nominal_fusion_twas_z_score[indices]))
					tissue_observed[tissue_iter] = 1.0
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
				raw_arr <- c(raw_arr, tissue_z[tissue_iter])
				observed_arr <- c(observed_arr, tissue_observed[tissue_iter])
			}

		}
	}
	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, tissue_causal_prob=count_arr, raw_statistic=raw_arr, observed=observed_arr)
	return(df)
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

extract_joint_susie_df <- function(trait_data_list, trait_names, tissue_names, column_name) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	raw_arr <- c()
	observed_arr <- c()
	gene_component_boolean_arr <- c()
	variant_overlap_prob_arr <- c()


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
			tissue_observed <- rep(0.0, num_tiss)
			component_has_gene_component <- 0.0
			tissue_variant_overlap_prob <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissue_observed[tissue_iter] = 1
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$tgfm_susie_pip[indices])) > .5)
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + max((component_df$variant_in_a_cs_prob[indices] > .5)*(component_df[column_name][,1])[indices])
					if (max((component_df[column_name][,1])[indices]) > 0.0) {
						component_has_gene_component = 1.0
					}
					best_col = which((component_df[column_name][,1])[indices]==max((component_df[column_name][,1])[indices]))
					tissue_variant_overlap_prob[tissue_iter] = component_df$variant_overlap_prob[indices][best_col]
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs((component_df[column_name][,1])[indices])) > .5)
				}
			}

			if (sum(tissues_un_normalized) ==0.0) {
				tissues_un_normalized_pseudo = tissues_un_normalized + pseudocount
			} else {
				tissues_un_normalized_pseudo = tissues_un_normalized
			}
			tissues_normalized =tissues_un_normalized_pseudo/sum(tissues_un_normalized_pseudo)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
				raw_arr <- c(raw_arr, tissues_un_normalized[tissue_iter])
				observed_arr <- c(observed_arr, tissue_observed[tissue_iter])
				gene_component_boolean_arr <- c(gene_component_boolean_arr, component_has_gene_component)
				variant_overlap_prob_arr <- c(variant_overlap_prob_arr, tissue_variant_overlap_prob[tissue_iter])
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, tissue_causal_prob=count_arr, raw_statistic=raw_arr, observed=observed_arr, gene_component_boolean=gene_component_boolean_arr, variant_overlap_prob=variant_overlap_prob_arr)
	return(df)
}

extract_susie_df <- function(trait_data_list, trait_names, tissue_names, column_name) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	raw_arr <- c()
	observed_arr <- c()
	gene_component_boolean_arr <- c()
	variant_overlap_prob_arr <- c()


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
			tissue_observed <- rep(0.0, num_tiss)
			component_has_gene_component <- 0.0
			tissue_variant_overlap_prob <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissue_observed[tissue_iter] = 1
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$tgfm_susie_pip[indices])) > .5)
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + max((component_df[column_name][,1])[indices])
					if (max((component_df[column_name][,1])[indices]) > 0.0) {
						component_has_gene_component = 1.0
					}
					best_col = which((component_df[column_name][,1])[indices]==max((component_df[column_name][,1])[indices]))[1]
					tissue_variant_overlap_prob[tissue_iter] = component_df$variant_overlap_prob[indices][best_col]
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs((component_df[column_name][,1])[indices])) > .5)
				}
			}

			if (sum(tissues_un_normalized) ==0.0) {
				tissues_un_normalized_pseudo = tissues_un_normalized + pseudocount
			} else {
				tissues_un_normalized_pseudo = tissues_un_normalized
			}
			tissues_normalized =tissues_un_normalized_pseudo/sum(tissues_un_normalized_pseudo)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
				raw_arr <- c(raw_arr, tissues_un_normalized[tissue_iter])
				observed_arr <- c(observed_arr, tissue_observed[tissue_iter])
				gene_component_boolean_arr <- c(gene_component_boolean_arr, component_has_gene_component)
				variant_overlap_prob_arr <- c(variant_overlap_prob_arr, tissue_variant_overlap_prob[tissue_iter])
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, tissue_causal_prob=count_arr, raw_statistic=raw_arr, observed=observed_arr, gene_component_boolean=gene_component_boolean_arr, variant_overlap_prob=variant_overlap_prob_arr)
	return(df)
}

extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_multivariate_susie_alpha_z <- function(trait_data_list, trait_names, tissue_names, column_name) {
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
			tissue_raw <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$tgfm_susie_pip[indices])) > .5)
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + max((component_df[column_name][,1])[indices])
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs((component_df[column_name][,1])[indices])) > .5)
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


extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_multivariate_alpha_z <- function(trait_data_list, trait_names, tissue_names) {
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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$tgfm_const_prior_twas_z_score[indices])) > 3.0)
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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$tgfm_tissue_prior_twas_z_score[indices])) >3.0)

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



extract_multivariate_z_robust_tissue_prior_df <- function(trait_data_list, trait_names, tissue_names, z_score_thresh) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	observed_arr <- c()
	raw_arr <- c()


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
			tissue_z <- rep(0.0, num_tiss)
			tissue_observed <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(abs(component_df$robust_tgfm_tissue_prior_twas_z_score[indices])) > z_score_thresh)
					tissue_z[tissue_iter] = max(abs(component_df$robust_tgfm_tissue_prior_twas_z_score[indices]))
					tissue_observed[tissue_iter] = 1.0
				}
			
			}

			tissues_un_normalized_pseudo = tissues_un_normalized + pseudocount
			tissues_normalized =tissues_un_normalized_pseudo/sum(tissues_un_normalized_pseudo)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
				raw_arr <- c(raw_arr, tissue_z[tissue_iter])
				observed_arr <- c(observed_arr, tissue_observed[tissue_iter])
			}

		}
	}
	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, tissue_causal_prob=count_arr, raw_statistic=raw_arr, observed=observed_arr)
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

extract_tgfm_sum_posterior_df <- function(trait_data_list, trait_names, tissue_names, column_name) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	observed_arr <- c()
	raw_arr <- c()


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
			tissue_observed <- rep(0, num_tiss)
			tissue_raw <- rep(0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum((component_df[column_name][,1])[indices])
					tissue_observed[tissue_iter] = 1
					tissue_raw[tissue_iter] = max(component_df$robust_tgfm_rss_regression_tissue_prior_log_likelihood[indices])
				}
			
			}

			tissues_normalized = tissues_un_normalized/sum(tissues_un_normalized)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
				observed_arr <- c(observed_arr, tissue_observed[tissue_iter])
				raw_arr <- c(raw_arr, tissue_raw[tissue_iter])
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, tissue_causal_prob=count_arr, raw_statistic=raw_arr, observed=observed_arr)

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

make_tissue_and_gene_set_variance_heatmaps_in_single_trait <- function(pseudotissue_gtex_rss_multivariate_twas_dir, trait_names, fusion_weights, robust, trait_name) {
	precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_False_robust_tissue_and_gene_specific_prior_precision_temp2.txt")

	df <- read.table(precision_file, header=TRUE, sep="\t")

	# PUt in inorganized data frame
	df$variance = 1.0/df$expected_precision

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	p <- ggplot(df, aes(x = tissue, y = gene_set, fill = variance)) +
  		geom_tile() +
  		theme(text = element_text(size=10), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		labs(fill="Expected Variance",title=paste0(trait_name))
  	return(p)

}

make_tissue_and_gene_set_variance_heatmaps <- function(pseudotissue_gtex_rss_multivariate_twas_dir, trait_names, fusion_weights, robust) {
	trait_arr <- c()
	gene_set_arr <- c()
	relative_variance_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_False_robust_tissue_and_gene_specific_prior_precision_temp2.txt")
		df <- read.table(precision_file, header=TRUE, sep="\t")
		df$variance = 1.0/df$expected_precision
		tissues <- as.character(sort(unique(df$tissue)))
		gene_sets <- as.character(sort(unique(df$gene_set)))
		#for (tissue_iter in 1:length(tissues)) {
			tissue <- tissues[2]
			df_tiss = df[df$tissue==tissue,]
			tiss_bgrd_var = df_tiss$variance[df_tiss$gene_set=="Universe"][1]
			for (gene_set_iter in 1:length(gene_sets)) {
				gene_set <- gene_sets[gene_set_iter]
				gene_set_var <- df_tiss$variance[df_tiss$gene_set==gene_set][1]
				ratio <- gene_set_var/tiss_bgrd_var
				relative_variance_arr <- c(relative_variance_arr, ratio)
				gene_set_arr <- c(gene_set_arr, gene_set)
				trait_arr <- c(trait_arr, trait_name)

			}
		#}
		#print(tissues)
		#print(gene_sets)

	}
	relative_variance_arr[relative_variance_arr > 15] = 15
	df_gene_set = data.frame(trait=factor(trait_arr, levels=trait_names), gene_set=factor(gene_set_arr, gene_sets), relative_variance=(relative_variance_arr))



	# PUt in inorganized data frame

	p_gene_set <- ggplot(df_gene_set, aes(x = trait, y = gene_set, fill = relative_variance)) +
  		geom_tile() +
  		theme(text = element_text(size=10), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		labs(fill="Expected Relative Variance") +
  		#scale_fill_gradientn(colours = c("red", "white", "blue"), values = c(0,1,3))
  		scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1.0)


	trait_arr <- c()
	tissue_arr <- c()
	variance_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		precision_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_cis_heritable_genes_count_genes_once_null_init_fusion_weights_False_robust_tissue_and_gene_specific_prior_precision_temp2.txt")
		df <- read.table(precision_file, header=TRUE, sep="\t")
		df$variance = 1.0/df$expected_precision
		tissues <- as.character(sort(unique(df$tissue)))
		gene_sets <- as.character(sort(unique(df$gene_set)))
		for (tissue_iter in 1:length(tissues)) {
			tissue <- tissues[tissue_iter]
			df_tiss = df[df$tissue==tissue,]
			tiss_bgrd_var = df_tiss$variance[df_tiss$gene_set=="Universe"][1]
			trait_arr <- c(trait_arr, trait_name)
			variance_arr <- c(variance_arr, tiss_bgrd_var)
			tissue_arr <- c(tissue_arr, tissue)

		}
	}
	df_tiss = data.frame(trait=factor(trait_arr, levels=trait_names), tissue=factor(tissue_arr, levels=tissues), variance=(variance_arr))

	df_tiss$tissue = str_replace_all(as.character(df_tiss$tissue), "-", "_")
	df_tiss$tissue <- recode(df_tiss$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	p_tiss <- ggplot(df_tiss, aes(x = tissue, y = trait, fill = variance)) +
  		geom_tile() +
  		theme(text = element_text(size=10), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		labs(fill="Expected Variance")

  	p <- plot_grid(p_gene_set, p_tiss, ncol=1, rel_heights=c(1, .7))
  	return(p)
}
tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z_v2 <- function(tgfm_robust_tissue_prior_df, susie_tissue_prior_df, multivariate_z_robust_tissue_prior_df, trait_name, single_causal_tissue) {
	tgfm_df <- tgfm_robust_tissue_prior_df[tgfm_robust_tissue_prior_df$tissue==single_causal_tissue,]
	susie_df <- susie_tissue_prior_df[susie_tissue_prior_df$tissue==single_causal_tissue,]
	z_df <- multivariate_z_robust_tissue_prior_df[multivariate_z_robust_tissue_prior_df$tissue==single_causal_tissue,]



	df <- data.frame(diff=tgfm_df$tissue_causal_prob-susie_df$tissue_causal_prob, twas_z=z_df$raw_statistic)

	p <- ggplot(df, aes(x=twas_z, y=diff)) +
  		geom_point(size=.7) +
  		figure_theme()

  	return(p)
}

tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z_v3 <- function(tgfm_robust_tissue_prior_df, susie_tissue_prior_df, multivariate_z_robust_tissue_prior_df, trait_name, single_causal_tissue) {
	tgfm_df <- tgfm_robust_tissue_prior_df[tgfm_robust_tissue_prior_df$tissue==single_causal_tissue,]
	susie_df <- susie_tissue_prior_df[susie_tissue_prior_df$tissue==single_causal_tissue,]
	z_df <- multivariate_z_robust_tissue_prior_df[multivariate_z_robust_tissue_prior_df$tissue==single_causal_tissue,]



	df <- data.frame(diff=tgfm_df$tissue_causal_prob-susie_df$tissue_causal_prob, variant_overlap_prob=susie_tissue_prior_df$variant_overlap_prob)

	p <- ggplot(df, aes(x=variant_overlap_prob, y=diff)) +
  		geom_point(size=.7) +
  		figure_theme()

  	return(p)
}




tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z <- function(tgfm_robust_tissue_prior_df, susie_tissue_prior_df, multivariate_z_robust_tissue_prior_df, trait_name, single_causal_tissue) {
	tgfm_df <- tgfm_robust_tissue_prior_df[tgfm_robust_tissue_prior_df$tissue==single_causal_tissue,]
	susie_df <- susie_tissue_prior_df[susie_tissue_prior_df$tissue==single_causal_tissue,]
	z_df <- multivariate_z_robust_tissue_prior_df[multivariate_z_robust_tissue_prior_df$tissue==single_causal_tissue,]



	df <- data.frame(tgfm_prob=tgfm_df$tissue_causal_prob, susie_prob=susie_df$tissue_causal_prob, twas_z=z_df$raw_statistic)

	p <- ggplot(df, aes(x=tgfm_prob, y=susie_prob, color=twas_z)) +
  		geom_point(position = position_jitter(w = 0.1, h = 0.1),size=.5) +
  		figure_theme()

  	return(p)
}

causal_gene_tissue_pairs_vs_variant_overlap_prob_boxplot <- function(df) {
	causal_indices <- df$tissue_causal_prob > .5
	not_causal_indices <- df$tissue_causal_prob <= .1


	variant_overlap_prob_arr = c(df$variant_overlap_prob[causal_indices], df$variant_overlap_prob[not_causal_indices])
	gene_class_arr <- c(rep("causal_gene", sum(causal_indices)), rep("not_causal_gene", sum(not_causal_indices)))

	variant_overlap_prob_arr[variant_overlap_prob_arr > .1] = .1

	df <- data.frame(variant_overlap_prob=variant_overlap_prob_arr, gene_class=gene_class_arr)


	p <- ggplot(df, aes(x=gene_class, y=variant_overlap_prob)) + 
  		geom_boxplot() +
  		figure_theme()

  	return(p)

}



make_number_of_susie_components_histogram <- function(num_components, model_name) {
	df = data.frame(num_components=num_components)
	p<-ggplot(df, aes(x=num_components)) + 
 	 geom_histogram(color="black", fill="white") +
 	 labs(x="Number of gene components/region", title=model_name) +
 	 scale_x_continuous(breaks=seq(0,10,1))+
 	 figure_theme()
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
trait_names_big <- as.character(trait_df$study_name)
print(trait_names_big)

for (itera in 1:length(trait_names_big)) {

	trait_names <- c(trait_names_big[itera])
	trait_name <- trait_names_big[itera]

	if (trait_name == "blood_WHITE_COUNT") {
		causal_tissue =c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "body_WHRadjBMIz") {
		causal_tissue = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum")
		single_causal_tissue = "Adipose_Subcutaneous"
	}
	if (trait_name == "biochemistry_Cholesterol") {
		causal_tissue = c("Liver")
		single_causal_tissue = "Liver"
	}
	if (trait_name == "bp_DIASTOLICadjMEDz") {
		causal_tissue = c("Artery_Heart", "Artery_Tibial", "Heart_Atrial_Appendage", "Heart_Left_Ventricle")
		single_causal_tissue = "Artery_Tibial"
	}
	if (trait_name == "blood_EOSINOPHIL_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "blood_RBC_DISTRIB_WIDTH") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "blood_RED_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "cov_EDU_COLLEGE") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
		single_causal_tissue = "Brain_Cerebellum"
	}
	if (trait_name == "mental_NEUROTICISM") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
		single_causal_tissue = "Brain_Cerebellum"
	}
	if (trait_name == "disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}


	fusion_weights="False"

	trait_to_tgfm_ch_data_hash = list()
	for (trait_iter in 1:length(trait_names)) {
		gene_version = "cis_heritable_genes"
		trait_name <- trait_names[trait_iter]
		trait_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_component_organized_with_and_without_genome_prior_tgfm_results.txt")
		df <- read.table(trait_file, header=TRUE)
		trait_to_tgfm_ch_data_hash[[trait_iter]] = df
	}



	number_susie_components_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_number_of_susie_components.txt")
	num_comp_df <- read.table(number_susie_components_file,header=TRUE)

	num_components_histo <- make_number_of_susie_components_histogram(num_comp_df$number_susie_components_tp, "tissue_prior_susie")
	num_components_histo2 <- make_number_of_susie_components_histogram(num_comp_df$number_susie_components_robust_tp, "robust_tissue_prior_susie_pip")
	merged_histo <- plot_grid(num_components_histo, num_components_histo2, ncol=1)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_number_of_susie_components.pdf")
	ggsave(merged_histo, file=output_file, width=7.2, height=6.0, units="in")



	susie_tissue_prior_df <- extract_susie_df(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "tgfm_tissue_prior_susie_pip")
	susie_robust_tissue_prior_df <- extract_susie_df(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_tissue_prior_susie_pip")
	#susie_joint_tissue_prior_df <- extract_joint_susie_df(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "tgfm_tissue_prior_susie_pip")
	tgfm_robust_tissue_prior_df <- extract_tgfm_sum_posterior_df(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob")
	multivariate_z_robust_tissue_prior_df_1 <- extract_multivariate_z_robust_tissue_prior_df(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, 1.0)
	multivariate_z_robust_tissue_prior_df_3 <- extract_multivariate_z_robust_tissue_prior_df(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, 3.0)
	univariate_z_df <- extract_univariate_z_df(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)





	prob_list <- list()
	prob_list[[1]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_univariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[2]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[3]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[4]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[5]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)


	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_fraction_causal_mediated_in_correct_tiss_twas_only.pdf")
	causal_tissue_prob_se_plot <- make_causal_tissue_prob(methods, prob_list, trait_name, length(causal_tissue)/38.0)
	ggsave(causal_tissue_prob_se_plot, file=output_file, width=7.2, height=5.0, units="in")




	#########################
	# Relationship between causal gene tissue pairs 
	#########################
	boxplot <- causal_gene_tissue_pairs_vs_variant_overlap_prob_boxplot(susie_tissue_prior_df)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_susie_causal_gene_tissue_pairs_vs_variant_overlap_prob_boxplot.pdf")
	ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")

	boxplot <- causal_gene_tissue_pairs_vs_variant_overlap_prob_boxplot(susie_robust_tissue_prior_df)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_rsusie_causal_gene_tissue_pairs_vs_variant_overlap_prob_boxplot.pdf")
	ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")

	#########################
	# Average causal prob bar plot
	#########################
	p_univariate <- make_average_causal_prob_bar_plot(univariate_z_df, "univariate const prior twas", causal_tissue, tissue_names)
	p_multivariate_1 <- make_average_causal_prob_bar_plot(multivariate_z_robust_tissue_prior_df_1, "mv tissue-prior r-TWAS (|zscore| > 1)", causal_tissue, tissue_names)
	p_multivariate_3 <- make_average_causal_prob_bar_plot(multivariate_z_robust_tissue_prior_df_3, "mv tissue-prior r-TWAS (|zscore| > 3)", causal_tissue, tissue_names)
	p_tgfm <- make_average_causal_prob_bar_plot(tgfm_robust_tissue_prior_df, "TGFM", causal_tissue, tissue_names)
	p_susie <- make_average_causal_prob_bar_plot(susie_tissue_prior_df, "tissue-prior susie-twas", causal_tissue, tissue_names)
	p_robust_susie <- make_average_causal_prob_bar_plot(susie_robust_tissue_prior_df, "tissue-prior r-susie-twas", causal_tissue, tissue_names)

	p_merged = plot_grid(p_univariate, p_tgfm, p_multivariate_1, p_multivariate_3, p_susie, p_robust_susie, ncol=2)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_fraction_causal_mediated.pdf")
	ggsave(p_merged, file=output_file, width=10.7, height=9.5, units="in")

	p_merged = plot_grid(p_tgfm, p_susie, ncol=1)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_fraction_causal_mediated_v2.pdf")
	ggsave(p_merged, file=output_file, width=7.2, height=7.0, units="in")

	#########################
	# Average causal prob bar plot for susie models for gwas components with at least one gene components
	#########################
	p_susie <- make_average_causal_prob_for_variant_components_with_nearby_gene_components_bar_plot(susie_tissue_prior_df, "tissue-prior susie-twas", causal_tissue, tissue_names)
	p_robust_susie <- make_average_causal_prob_for_variant_components_with_nearby_gene_components_bar_plot(susie_robust_tissue_prior_df, "tissue-prior r-susie-twas", causal_tissue, tissue_names)
	p_merged = plot_grid(p_susie, p_robust_susie, ncol=1)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_fraction_causal_mediated_among_variant_components_with_gene_components.pdf")
	ggsave(p_merged, file=output_file, width=7.2, height=7.0, units="in")


	#########################
	# TGFM causal prob vs Susie causal prob colored by twas Z-score
	#########################
	pp <- tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z(tgfm_robust_tissue_prior_df, susie_tissue_prior_df, multivariate_z_robust_tissue_prior_df_1, trait_name, single_causal_tissue)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_tgfm_vs_susie_scatter_colored_by_z.pdf")
	ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")

	pp <- tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z(tgfm_robust_tissue_prior_df, susie_robust_tissue_prior_df, multivariate_z_robust_tissue_prior_df_1, trait_name, single_causal_tissue)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_tgfm_vs_rsusie_scatter_colored_by_z.pdf")
	ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")

	pp <- tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z_v2(tgfm_robust_tissue_prior_df, susie_tissue_prior_df, univariate_z_df, trait_name, single_causal_tissue)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_tgfm_vs_susie_scatter_colored_by_z_v2.pdf")
	ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")

	pp <- tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z_v2(tgfm_robust_tissue_prior_df, susie_robust_tissue_prior_df, univariate_z_df, trait_name, single_causal_tissue)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_tgfm_vs_rsusie_scatter_colored_by_z_v2.pdf")
	ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")

	pp <- tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z_v3(tgfm_robust_tissue_prior_df, susie_tissue_prior_df, univariate_z_df, trait_name, single_causal_tissue)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_tgfm_vs_susie_scatter_colored_by_z_v3.pdf")
	ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")

	pp <- tgfm_causal_prob_vs_susie_causal_prob_scatter_colored_by_twas_z_v3(tgfm_robust_tissue_prior_df, susie_robust_tissue_prior_df, univariate_z_df, trait_name, single_causal_tissue)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "v2_", trait_name ,"_tgfm_vs_rsusie_scatter_colored_by_z_v3.pdf")
	ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")

	print("DONE")



}



if (FALSE) {
proportion_trait_components_mediated_tgfm_univariate_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_univariate_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('one')
proportion_trait_components_mediated_tgfm_multivariate_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_multivariate_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('two')
proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_multivariate_susie_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('two')
proportion_trait_components_mediated_tgfm_rss_regression_const_prior_sum_posterior_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "tgfm_rss_regression_const_prior_posterior_prob")
print('two')
proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_tissue_prior_multivariate_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('three')
proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z_list = extract_proportion_of_trait_components_mediated_in_each_tissue_robust_tgfm_tissue_prior_multivariate_alpha_z(trait_to_tgfm_ch_data_hash, trait_names, tissue_names)
print('four')
proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "tgfm_rss_regression_tissue_prior_posterior_prob")
print('four')
proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list = extract_proportion_of_trait_components_mediated_in_each_tissue_tgfm_sum_posterior(trait_to_tgfm_ch_data_hash, trait_names, tissue_names, "robust_tgfm_rss_regression_tissue_prior_posterior_prob")
print('five')




proportion_trait_components_mediated_tgfm_univariate_alpha_z = proportion_trait_components_mediated_tgfm_univariate_alpha_z_list[[1]]
proportion_trait_components_mediated_tgfm_multivariate_alpha_z = proportion_trait_components_mediated_tgfm_multivariate_alpha_z_list[[1]]
proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z = proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z_list[[1]]
proportion_trait_components_mediated_tgfm_rss_regression_const_prior_sum_posterior = proportion_trait_components_mediated_tgfm_rss_regression_const_prior_sum_posterior_list[[1]]
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
	if (trait_name == "disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}

	p_univariate <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_univariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_univariate_alpha_z$trait)==trait_name,], "univariate const prior twas", causal_tissue)
	p_multivariate_cp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_multivariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_multivariate_alpha_z$trait)==trait_name,], "mv const prior TWAS", causal_tissue)
	p_susie_multivariate_cp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z[as.character(proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z$trait)==trait_name,], "mv const prior SuSiE-TWAS", causal_tissue)
	p_tgfm_cp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_rss_regression_const_prior_sum_posterior[as.character(proportion_trait_components_mediated_tgfm_rss_regression_const_prior_sum_posterior$trait)==trait_name,], "const prior TGFM", causal_tissue)
	p_multivariate_tp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z$trait)==trait_name,], "mv tissue-prior TWAS", causal_tissue)
	p_robust_multivariate_tp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z$trait)==trait_name,], "mv tissue-prior r-TWAS", causal_tissue)
	p_tgfm <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior[as.character(proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior$trait)==trait_name,], "tissue-prior TGFM", causal_tissue)
	p_robust_tgfm <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior[as.character(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior$trait)==trait_name,], "tissue-prior r-TGFM", causal_tissue)


	p_merged = plot_grid(p_univariate, p_multivariate_cp, p_susie_multivariate_cp, p_tgfm_cp, p_multivariate_tp, p_robust_multivariate_tp, p_tgfm, p_robust_tgfm, ncol=2)

	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_fraction_causal_mediated.pdf")
	ggsave(p_merged, file=output_file, width=10.7, height=9.5, units="in")
}



for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	methods <- c("univariate const prior twas", "mv const prior TWAS", "mv const prior SuSiE-TWAS", "const prior TGFM", "mv tissue-prior TWAS", "mv tissue-prior r-TWAS", "tissue-prior TGFM", "tissue-prior r-TGFM")

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
	if (trait_name == "disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}


	prob_list <- list()
	prob_list[[1]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_univariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[2]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[3]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[4]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_rss_regression_const_prior_sum_posterior_list[[2]], trait_name, causal_tissue)
	prob_list[[5]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[6]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[7]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_rss_regression_tissue_prior_sum_posterior_list[[2]], trait_name, causal_tissue)
	prob_list[[8]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_robust_tgfm_rss_regression_tissue_prior_sum_posterior_list[[2]], trait_name, causal_tissue)


	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_fraction_causal_mediated_in_correct_tiss.pdf")
	causal_tissue_prob_se_plot <- make_causal_tissue_prob(methods, prob_list, trait_name, length(causal_tissue)/38.0)
	ggsave(causal_tissue_prob_se_plot, file=output_file, width=7.2, height=5.0, units="in")
}





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
	if (trait_name == "disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}

	p_univariate <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_univariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_univariate_alpha_z$trait)==trait_name,], "univariate const prior twas", causal_tissue)
	p_multivariate_cp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_multivariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_multivariate_alpha_z$trait)==trait_name,], "mv const prior TWAS", causal_tissue)
	p_susie_multivariate_cp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z[as.character(proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z$trait)==trait_name,], "mv const prior SuSiE-TWAS", causal_tissue)
	p_multivariate_tp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z$trait)==trait_name,], "mv tissue-prior TWAS", causal_tissue)
	p_robust_multivariate_tp <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z[as.character(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z$trait)==trait_name,], "mv tissue-prior r-TWAS", causal_tissue)


	p_merged = plot_grid(p_univariate, p_multivariate_cp, p_susie_multivariate_cp, p_multivariate_tp, p_robust_multivariate_tp, ncol=2)

	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_fraction_causal_mediated_twas_only.pdf")
	ggsave(p_merged, file=output_file, width=10.7, height=9.5, units="in")
}



for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	methods <- c("univariate const prior twas", "mv const prior TWAS", "mv const prior SuSiE-TWAS", "mv tissue-prior TWAS", "mv tissue-prior r-TWAS")

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
	if (trait_name == "disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
	}


	prob_list <- list()
	prob_list[[1]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_univariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[2]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[3]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_multivariate_susie_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[4]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)
	prob_list[[5]] = get_probs_in_causal_tissues_cross_traits(proportion_trait_components_mediated_robust_tgfm_tissue_prior_multivariate_alpha_z_list[[2]], trait_name, causal_tissue)


	output_file <- paste0(rss_multivariate_twas_visualization_dir, trait_name,"_fraction_causal_mediated_in_correct_tiss_twas_only.pdf")
	causal_tissue_prob_se_plot <- make_causal_tissue_prob(methods, prob_list, trait_name, length(causal_tissue)/38.0)
	ggsave(causal_tissue_prob_se_plot, file=output_file, width=7.2, height=5.0, units="in")
}








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

