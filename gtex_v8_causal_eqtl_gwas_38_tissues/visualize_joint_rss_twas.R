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

make_average_causal_prob_bar_plot_across_tissues <- function(df_causal_prob, model_name, causal_tissue, tissue_names) {
	tissue_arr <- c()
	average_causal_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		if (model_name == "TGFM") {
			avg_prob <- sum(df_causal_prob$raw_statistic[as.character(df_causal_prob$tissue)==tissue_name])/sum(as.character(df_causal_prob$tissue)==tissue_name)
		}
		if (model_name == "Coloc") {
			avg_prob <- sum(df_causal_prob$tissue_causal_med_prob[as.character(df_causal_prob$tissue)==tissue_name])/sum(as.character(df_causal_prob$tissue)==tissue_name)
		}
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
  			labs(y="Fraction of\n trait components", x="",title=model_name) +
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


extract_coloc_df <- function(df, trait_name, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	raw_arr <- c()
	observed_arr <- c()
	tissue_causal_med_prob_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005

	#for (trait_iter in 1:length(trait_names)) {
		#trait_name <- as.character(trait_names[trait_iter])

		#df <- trait_data_list[[trait_iter]]

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
			tissue_causal_med_prob <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + (max(abs(component_df$coloc_pph4[indices])))
					tissue_observed[tissue_iter] = 1.0
				}
			
			}
			if (sum(tissues_un_normalized) == 0.0) {
				tissues_un_normalized = tissues_un_normalized + pseudocount
			} else {
				tissue_causal_med_prob = max(tissues_un_normalized)*(tissues_un_normalized/sum(tissues_un_normalized))
			}
			tissues_normalized =tissues_un_normalized/sum(tissues_un_normalized)


			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
				raw_arr <- c(raw_arr, tissues_un_normalized[tissue_iter])
				observed_arr <- c(observed_arr, tissue_observed[tissue_iter])
				tissue_causal_med_prob_arr <- c(tissue_causal_med_prob_arr, tissue_causal_med_prob[tissue_iter])
			}

		}
	#}
	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, tissue_causal_prob=count_arr, tissue_causal_med_prob=tissue_causal_med_prob_arr, raw_statistic=raw_arr, observed=observed_arr)
	return(df)
}

extract_univariate_z_df <- function(df, trait_name, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	raw_arr <- c()
	observed_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005

	#for (trait_iter in 1:length(trait_names)) {
		#trait_name <- as.character(trait_names[trait_iter])

		#df <- trait_data_list[[trait_iter]]

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
	#}
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

extract_joint_susie_non_med_prob_df <- function(df, trait_name, tissue_names, column_name) {
	trait_arr <- c()
	comp_arr <- c()
	non_med_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005



		
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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum((component_df[column_name][,1][indices]))
				}
			}

			trait_arr <- c(trait_arr, trait_name)
			comp_arr <- c(comp_arr, trait_component)
			non_med_arr <- c(non_med_arr, 1.0-sum(tissues_un_normalized))
		}
	
	df <- data.frame(trait=trait_arr, trait_component=comp_arr, non_mediated_probability=non_med_arr)
	return(df)
}


extract_joint_susie_df <- function(df, trait_name, tissue_names, column_name) {
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

	#for (trait_iter in 1:length(trait_names)) {
		#trait_name <- as.character(trait_names[trait_iter])

		#df <- trait_data_list[[trait_iter]]

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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + sum((component_df[column_name][,1][indices]))
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + max((component_df$variant_in_a_cs_prob[indices] > .5)*(component_df[column_name][,1])[indices])
					#if (max((component_df[column_name][,1])[indices]) > 0.0) {
					#	component_has_gene_component = 1.0
					#}
					#best_col = which((component_df[column_name][,1])[indices]==max((component_df[column_name][,1])[indices]))
					#tissue_variant_overlap_prob[tissue_iter] = component_df$variant_overlap_prob[indices][best_col]
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
				#gene_component_boolean_arr <- c(gene_component_boolean_arr, component_has_gene_component)
			}

		}
	#}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, tissue_causal_prob=count_arr, raw_statistic=raw_arr, observed=observed_arr)
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



extract_multivariate_z_robust_tissue_prior_df <- function(df, trait_name, tissue_names, z_score_thresh) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	observed_arr <- c()
	raw_arr <- c()


	num_tiss <- length(tissue_names)

	pseudocount = .005

	#for (trait_iter in 1:length(trait_names)) {
		#trait_name <- as.character(trait_names[trait_iter])

		#df <- trait_data_list[[trait_iter]]

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
	#}
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

extract_tgfm_sum_posterior_df <- function(df, trait_name, tissue_names, column_name) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	observed_arr <- c()
	raw_arr <- c()


	num_tiss <- length(tissue_names)



	#for (trait_iter in 1:length(trait_names)) {
		#trait_name <- as.character(trait_names[trait_iter])

		#df <- trait_data_list[[trait_iter]]
		
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
	#}

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
	return(p)
}


make_distribution_of_mediated_probs_in_each_tissue <- function(df, trait_name, tissue_names) {
	tissue_names <- as.character(tissue_names)

	tissue_arr <- c()
	posterior_prob_arr <- c()
	causal_tissue_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]

		probs <- df$raw_statistic[as.character(df$tissue)==tissue_name]

		posterior_prob_arr <- c(posterior_prob_arr, probs)
		tissue_arr <- c(tissue_arr, rep(tissue_name, length(probs)))
	}

	df <- data.frame(tissue=tissue_arr, posterior_prob=posterior_prob_arr)
	p <- ggplot(df, aes(x=posterior_prob)) + geom_histogram() + geom_histogram() + figure_theme() +
		facet_wrap( ~ tissue, ncol = 4) + 
		labs(x="Expression-mediated probability / component")
	return(p)
}


make_distribution_of_coloc_probs_in_each_tissue <- function(df, trait_name, tissue_names) {
	tissue_names <- as.character(tissue_names)

	tissue_arr <- c()
	posterior_prob_arr <- c()
	causal_tissue_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]

		probs <- df$raw_statistic[as.character(df$tissue)==tissue_name]


		posterior_prob_arr <- c(posterior_prob_arr, probs)
		tissue_arr <- c(tissue_arr, rep(tissue_name, length(probs)))
	}

	df <- data.frame(tissue=tissue_arr, posterior_prob=posterior_prob_arr)
	p <- ggplot(df, aes(x=posterior_prob)) + geom_histogram() + geom_histogram() + figure_theme() +
		facet_wrap( ~ tissue, ncol = 4) + 
		labs(x="Max coloc pph4 / component")
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

summarize_mediated_components_for_joint_susie_model_for_pres <- function(susie_joint_tissue_prior_df, tissue_names, trait_name, gtex_colors_df) {
	# Standard summary
	tissue_arr <- c()
	fraction_arr <- c()
	num_components <- length(unique(susie_joint_tissue_prior_df$trait_component))

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		fraction <- sum(susie_joint_tissue_prior_df$raw_statistic[as.character(susie_joint_tissue_prior_df$tissue)==tissue_name])/num_components

		tissue_arr <- c(tissue_arr, paste0(tissue_name))
		fraction_arr <- c(fraction_arr, fraction)
	}
	
	non_med_fraction <- 1.0 - (sum(susie_joint_tissue_prior_df$raw_statistic)/num_components)
	tissue_arr <- c(tissue_arr, "Non-mediated")
	fraction_arr <- c(fraction_arr, non_med_fraction)
	df <- data.frame(tissue=tissue_arr, fraction=fraction_arr)

	# Summarize by top 3 tissues
	ordered_tissues <- tissue_arr[order(-fraction_arr)]


	tissue_arr2 <- c()
	fraction_arr2 <- c()
	trait_arr <- c()
	colors <- c()
	num_tissues <- 5
	for (indexer in 2:num_tissues) {
		if (ordered_tissues[indexer] == "Non-mediated") {
			colors <- c(colors, "#999999")
		} else {
			tt <- ordered_tissues[indexer]
			if (startsWith(tt, "Brain")) {
				colors <- c(colors, "#EEEE00")
			} else {
				index = which(gtex_colors_df$tissue_site_detail_id==tt)
				colors <- c(colors, paste0('#', as.character(gtex_colors_df$tissue_color_hex)[index]))
			}
		}

		tissue_arr2 <- c(tissue_arr2, paste0(ordered_tissues[indexer], "-mediated"))
		fraction_arr2 <- c(fraction_arr2, as.numeric(df$fraction[df$tissue == ordered_tissues[indexer]]))
		trait_arr <- c(trait_arr, trait_name)
	}

	tissue_arr2 <- c(tissue_arr2, "Other-expression-mediated")
	fraction_arr2 <- c(fraction_arr2, 1.0 - sum(fraction_arr2) - as.numeric(df$fraction[df$tissue == ordered_tissues[1]]))
	trait_arr <- c(trait_arr, trait_name)
	colors <- c(colors, "darkblue")

	for (indexer in 1:1) {
		if (ordered_tissues[indexer] == "Non-mediated") {
			colors <- c(colors, "#CCCCCC")
			namer <- "Non-mediated"
		} else {
			tt <- ordered_tissues[indexer]
			if (startsWith(tt, "Brain")) {
				colors <- c(colors, "#EEEE00")
			} else {
				index = which(gtex_colors_df$tissue_site_detail_id==tt)
				colors <- c(colors, paste0('#', as.character(gtex_colors_df$tissue_color_hex)[index]))
			}
			namer <- paste0(ordered_tissues[indexer], "-mediated")
		}

		tissue_arr2 <- c(tissue_arr2, namer)
		fraction_arr2 <- c(fraction_arr2, as.numeric(df$fraction[df$tissue == ordered_tissues[indexer]]))
		trait_arr <- c(trait_arr, trait_name)
	}



	# Compact df
	df2 <- data.frame(tissue=factor(tissue_arr2, levels=rev(c(tissue_arr2))), fraction=fraction_arr2, trait=trait_arr, color=colors)

	df2 <- df2[df2$tissue!="Non-mediated",]
	#colors <- as.character(df2$color)

	p <- ggplot(df2, aes(fill=tissue, y=trait, x=fraction)) + 
    	geom_bar(position="fill", stat="identity") +
    	figure_theme() +
    	labs(x="Fraction of expression-mediated disease components", title=trait_name, y="", fill="") +
    	scale_fill_manual(values=rev(as.character(df2$color))) +
    	theme(legend.position="bottom") +
    	theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank())  +
    	guides(fill = guide_legend(reverse=TRUE, nrow=2, byrow=TRUE)) +
    	theme(legend.key.size = unit(.5, 'cm'), legend.title = element_text(size=11),legend.text = element_text(size=11)) +
    	theme(plot.title = element_text(hjust = 0.5))
   	return(p)


}

summarize_mediated_components_for_joint_susie_model <- function(susie_joint_tissue_prior_df, tissue_names, trait_name, gtex_colors_df) {
	# Standard summary
	tissue_arr <- c()
	fraction_arr <- c()
	num_components <- length(unique(susie_joint_tissue_prior_df$trait_component))

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		fraction <- sum(susie_joint_tissue_prior_df$raw_statistic[as.character(susie_joint_tissue_prior_df$tissue)==tissue_name])/num_components

		tissue_arr <- c(tissue_arr, paste0(tissue_name))
		fraction_arr <- c(fraction_arr, fraction)
	}
	
	non_med_fraction <- 1.0 - (sum(susie_joint_tissue_prior_df$raw_statistic)/num_components)
	tissue_arr <- c(tissue_arr, "Non-mediated")
	fraction_arr <- c(fraction_arr, non_med_fraction)
	df <- data.frame(tissue=tissue_arr, fraction=fraction_arr)

	# Summarize by top 3 tissues
	ordered_tissues <- tissue_arr[order(-fraction_arr)]


	tissue_arr2 <- c()
	fraction_arr2 <- c()
	trait_arr <- c()
	colors <- c()
	num_tissues <- 5
	for (indexer in 2:num_tissues) {
		if (ordered_tissues[indexer] == "Non-mediated") {
			colors <- c(colors, "#999999")
		} else {
			tt <- ordered_tissues[indexer]
			if (startsWith(tt, "Brain")) {
				colors <- c(colors, "#EEEE00")
			} else {
				index = which(gtex_colors_df$tissue_site_detail_id==tt)
				colors <- c(colors, paste0('#', as.character(gtex_colors_df$tissue_color_hex)[index]))
			}
		}

		tissue_arr2 <- c(tissue_arr2, paste0(ordered_tissues[indexer], "-mediated"))
		fraction_arr2 <- c(fraction_arr2, as.numeric(df$fraction[df$tissue == ordered_tissues[indexer]]))
		trait_arr <- c(trait_arr, trait_name)
	}

	tissue_arr2 <- c(tissue_arr2, "Other-expression-mediated")
	fraction_arr2 <- c(fraction_arr2, 1.0 - sum(fraction_arr2) - as.numeric(df$fraction[df$tissue == ordered_tissues[1]]))
	trait_arr <- c(trait_arr, trait_name)
	colors <- c(colors, "darkblue")

	for (indexer in 1:1) {
		if (ordered_tissues[indexer] == "Non-mediated") {
			colors <- c(colors, "#CCCCCC")
			namer <- "Non-mediated"
		} else {
			tt <- ordered_tissues[indexer]
			if (startsWith(tt, "Brain")) {
				colors <- c(colors, "#EEEE00")
			} else {
				index = which(gtex_colors_df$tissue_site_detail_id==tt)
				colors <- c(colors, paste0('#', as.character(gtex_colors_df$tissue_color_hex)[index]))
			}
			namer <- paste0(ordered_tissues[indexer], "-mediated")
		}

		tissue_arr2 <- c(tissue_arr2, namer)
		fraction_arr2 <- c(fraction_arr2, as.numeric(df$fraction[df$tissue == ordered_tissues[indexer]]))
		trait_arr <- c(trait_arr, trait_name)
	}



	# Compact df
	df2 <- data.frame(tissue=factor(tissue_arr2, levels=rev(c(tissue_arr2))), fraction=fraction_arr2, trait=trait_arr, color=colors)

	df2 <- df2[df2$tissue!="Non-mediated",]
	#colors <- as.character(df2$color)

	p <- ggplot(df2, aes(fill=tissue, y=trait, x=fraction)) + 
    	geom_bar(position="fill", stat="identity") +
    	figure_theme() +
    	labs(x="Fraction of expression-mediated disease components", title=trait_name, y="", fill="") +
    	scale_fill_manual(values=rev(as.character(df2$color))) +
    	theme(legend.position="bottom") +
    	theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank())  +
    	guides(fill = guide_legend(reverse=TRUE, nrow=3, byrow=TRUE)) +
    	theme(legend.key.size = unit(.2, 'cm'), legend.title = element_text(size=9),legend.text = element_text(size=9)) +
    	theme(plot.title = element_text(hjust = 0.5))
   	return(p)
}

summarize_mediated_and_non_mediated_components_for_joint_susie_model <- function(susie_joint_tissue_prior_df, tissue_names, trait_name, gtex_colors_df) {
	# Standard summary
	tissue_arr <- c()
	fraction_arr <- c()
	num_components <- length(unique(susie_joint_tissue_prior_df$trait_component))

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		fraction <- sum(susie_joint_tissue_prior_df$raw_statistic[as.character(susie_joint_tissue_prior_df$tissue)==tissue_name])/num_components

		tissue_arr <- c(tissue_arr, paste0(tissue_name))
		fraction_arr <- c(fraction_arr, fraction)
	}
	
	non_med_fraction <- 1.0 - (sum(susie_joint_tissue_prior_df$raw_statistic)/num_components)
	tissue_arr <- c(tissue_arr, "Non-mediated")
	fraction_arr <- c(fraction_arr, non_med_fraction)
	df <- data.frame(tissue=tissue_arr, fraction=fraction_arr)

	# Summarize by top 3 tissues
	ordered_tissues <- tissue_arr[order(-fraction_arr)]


	tissue_arr2 <- c()
	fraction_arr2 <- c()
	trait_arr <- c()
	colors <- c()
	num_tissues <- 3
	for (indexer in 2:num_tissues) {
		if (ordered_tissues[indexer] == "Non-mediated") {
			colors <- c(colors, "#999999")
			namer <- "Non-mediated"
		} else {
			tt <- ordered_tissues[indexer]
			if (startsWith(tt, "Brain")) {
				colors <- c(colors, "#EEEE00")
			} else {
				index = which(gtex_colors_df$tissue_site_detail_id==tt)
				colors <- c(colors, paste0('#', as.character(gtex_colors_df$tissue_color_hex)[index]))
			}
			namer <- paste0(ordered_tissues[indexer], "-mediated")
		}

		tissue_arr2 <- c(tissue_arr2, namer)
		fraction_arr2 <- c(fraction_arr2, as.numeric(df$fraction[df$tissue == ordered_tissues[indexer]]))
		trait_arr <- c(trait_arr, trait_name)
	}

	tissue_arr2 <- c(tissue_arr2, "Other-expression-mediated")
	fraction_arr2 <- c(fraction_arr2, 1.0 - sum(fraction_arr2) - as.numeric(df$fraction[df$tissue == ordered_tissues[1]]))
	trait_arr <- c(trait_arr, trait_name)
	colors <- c(colors, "darkblue")

	for (indexer in 1:1) {
		if (ordered_tissues[indexer] == "Non-mediated") {
			colors <- c(colors, "#CCCCCC")
			namer <- "Non-mediated"
		} else {
			tt <- ordered_tissues[indexer]
			if (startsWith(tt, "Brain")) {
				colors <- c(colors, "#EEEE00")
			} else {
				index = which(gtex_colors_df$tissue_site_detail_id==tt)
				colors <- c(colors, paste0('#', as.character(gtex_colors_df$tissue_color_hex)[index]))
			}
			namer <- paste0(ordered_tissues[indexer], "-mediated")
		}

		tissue_arr2 <- c(tissue_arr2, namer)
		fraction_arr2 <- c(fraction_arr2, as.numeric(df$fraction[df$tissue == ordered_tissues[indexer]]))
		trait_arr <- c(trait_arr, trait_name)
	}



	# Compact df
	df2 <- data.frame(tissue=factor(tissue_arr2, levels=rev(c(tissue_arr2))), fraction=fraction_arr2, trait=trait_arr, color=colors)
	#print(df2)

	p <- ggplot(df2, aes(fill=tissue, y=trait, x=fraction)) + 
    	geom_bar(position="fill", stat="identity") +
    	figure_theme() +
    	labs(x="Fraction of disease components", title=trait_name, y="", fill="") +
    	scale_fill_manual(values=rev(colors)) +
    	theme(legend.position="bottom") +
    	theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank())  +
    	guides(fill = guide_legend(reverse=TRUE, nrow=3, byrow=TRUE)) +
    	theme(legend.key.size = unit(.2, 'cm'), legend.title = element_text(size=9),legend.text = element_text(size=9)) +
    	theme(plot.title = element_text(hjust = 0.5))
   	return(p)
}


make_mediated_prob_se_barplot <- function(df) {
	p <- ggplot(df) +
    		geom_bar( aes(x=trait, y=average_mediated_probability), stat="identity", fill="skyblue", alpha=0.7) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    		labs(y="Average expression-mediated probability", x="") +
    		geom_errorbar( aes(x=trait, ymin=average_mediated_probability-(1.96*average_mediated_probability_se), ymax=average_mediated_probability+(1.96*average_mediated_probability_se)), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    		figure_theme()
    return(p)
}

make_histogram_showing_distribution_of_mediated_probabilities_across_components <- function(df, trait_name) {
	df$mediated_probability = 1.0 - df$non_mediated_probability

	p <- ggplot(df, aes(x=mediated_probability))+
	  geom_histogram(color="darkblue", fill="lightblue") +
	  figure_theme() +
	  labs(x="Expression-mediated probability", y="Number of components", title=trait_name) +
	  theme(plot.title = element_text(hjust = 0.5))

	return(p)

}


get_causal_tissues_corresponding_to_trait <- function(trait_name) {
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
		single_causal_tissue = "Cells_EBV_transformed_lymphocytes"
	}

	ret <- list()
	ret[[1]] = causal_tissue
	ret[[2]] = single_causal_tissue

	return(ret)
}

trait_file <- args[1]
gtex_pseudotissue_file <- args[2]
pseudotissue_gtex_rss_multivariate_twas_dir <- args[3]
rss_multivariate_twas_visualization_dir <- args[4]
gtex_tissue_colors_file <- args[5]

if (FALSE) {
# Load in gtex tissues
tissue_df <- read.table(gtex_pseudotissue_file,header=TRUE)
tissue_names <- as.character(tissue_df$pseudotissue_name)
tissue_names <- str_replace_all(tissue_names, "-", "_")

# Load in gtex tissue colors 
gtex_colors_df <- read.table(gtex_tissue_colors_file, header=TRUE, sep="\t")
gtex_colors_df$tissue_site_detail_id = as.character(gtex_colors_df$tissue_site_detail_id)
gtex_colors_df$tissue_site_detail_id[23] = "Cells_Cultured_fibroblasts"


# trait names
trait_df <- read.table(trait_file, header=TRUE)
trait_df = trait_df[as.character(trait_df$study_name) != "biochemistry_Triglycerides", ]
trait_df = trait_df[as.character(trait_df$study_name) != "mental_NEUROTICISM", ]
trait_df = trait_df[as.character(trait_df$study_name) != "disease_ALLERGY_ECZEMA_DIAGNOSED", ]

trait_names_big <- as.character(trait_df$study_name)

# Model parameters
fusion_weights="False"
gene_version = "cis_heritable_genes"


# Keep track of average mediated probability and se
avg_med_prob_arr <- c()
se_med_prob_arr <- c()


# Create list to keep track of stacked barplot summary plots across traits
stacked_barplot_joint_susie_list <- list()
stacked_barplot_joint_susie_mediated_list <- list()

# Loop through traits
for (itera in 1:length(trait_names_big)) {

	# Name of trait defining current loop
	trait_name <- trait_names_big[itera]

	# Get causal tissue and single_causal tissue corresponding to this trait
	causal_tissue_list <- get_causal_tissues_corresponding_to_trait(trait_name)
	causal_tissue <- causal_tissue_list[[1]]
	single_causal_tissue <- causal_tissue_list[[2]]

	# Load in organized TGFM results
	trait_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_component_organized_joint_tgfm_results.txt")
	trait_df <- read.table(trait_file, header=TRUE)
	

	# Calculate non mediated probability across all components for this trait
	susie_joint_non_med_prob_df <- extract_joint_susie_non_med_prob_df(trait_df, trait_name, tissue_names, "robust_tgfm_tissue_prior_joint_susie_prob")
	susie_joint_tissue_prior_df <- extract_joint_susie_df(trait_df, trait_name, tissue_names, "robust_tgfm_tissue_prior_joint_susie_prob")
	non_med_prob = 1.0 - (sum(susie_joint_tissue_prior_df$raw_statistic)/length(unique(susie_joint_tissue_prior_df$trait_component)))
	med_prob = 1.0 - non_med_prob


	# Calculate probabilities generated by other methods
	coloc_df <- extract_coloc_df(trait_df, trait_name, tissue_names)

	# Histogram showing distribution of mediated probabilities across components
	histo <- make_histogram_showing_distribution_of_mediated_probabilities_across_components(susie_joint_non_med_prob_df, trait_name)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "mediated_probability_histogram_", trait_name, ".pdf")
	ggsave(histo, file=output_file, width=7.2, height=4.0, units="in")



	# Summarize mediated and non-mediated probabilities across components
	stacked_barplot_summary <- summarize_mediated_and_non_mediated_components_for_joint_susie_model(susie_joint_tissue_prior_df, tissue_names, trait_name, gtex_colors_df) 
	stacked_barplot_joint_susie_list[[itera]] = stacked_barplot_summary
	#output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_susie_stacked_barplot_summary_", trait_name, ".pdf")
	#ggsave(stacked_barplot_summary, file=output_file, width=7.2, height=4.0, units="in")


	# Summarize mediated and probabilities across components
	stacked_barplot_mediated_summary <- summarize_mediated_components_for_joint_susie_model(susie_joint_tissue_prior_df, tissue_names, trait_name, gtex_colors_df) 
	stacked_barplot_joint_susie_mediated_list[[itera]] = stacked_barplot_mediated_summary
	#output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_susie_mediated_components_stacked_barplot_summary_", trait_name, ".pdf")
	#ggsave(stacked_barplot_mediated_summary, file=output_file, width=7.2, height=4.0, units="in")


	# Keep track of average mediated probability (and its standard errror) for each trait
	avg_med_prob_arr <- c(avg_med_prob_arr, mean(1.0 - susie_joint_non_med_prob_df$non_mediated_probability))
	se_med_prob_arr <- c(se_med_prob_arr, std_mean(1.0 - susie_joint_non_med_prob_df$non_mediated_probability))



	# Make distribution of mediated probabilities by each tissue
	med_prob_histo <- make_distribution_of_mediated_probs_in_each_tissue(susie_joint_tissue_prior_df, trait_name, tissue_names)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_susie_mediated_probability_per_tissue_", trait_name, ".pdf")
	ggsave(med_prob_histo, file=output_file, width=9.2, height=11.0, units="in")

	# Make distribution of coloc probabilities by each tissue
	coloc_prob_histo <- make_distribution_of_coloc_probs_in_each_tissue(coloc_df, trait_name, tissue_names)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "coloc_probability_per_tissue_", trait_name, ".pdf")
	ggsave(coloc_prob_histo, file=output_file, width=9.2, height=11.0, units="in")


	# Show mean number of mediated components in each tissue
	tgfm_med_comp_per_tissue_bar <- make_average_causal_prob_bar_plot_across_tissues(susie_joint_tissue_prior_df, "TGFM", causal_tissue, tissue_names)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_susie_mediated_component_per_tissue_bar_plot_", trait_name, ".pdf")
	ggsave(tgfm_med_comp_per_tissue_bar, file=output_file, width=7.2, height=4.0, units="in")
	# Show mean number of mediated components in each tissue
	coloc_med_comp_per_tissue_bar <- make_average_causal_prob_bar_plot_across_tissues(coloc_df, "Coloc", causal_tissue, tissue_names)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "coloc_mediated_component_per_tissue_bar_plot_", trait_name, ".pdf")
	ggsave(coloc_med_comp_per_tissue_bar, file=output_file, width=7.2, height=4.0, units="in")
	####
	#JOINT VERSION
	####
	joint_med_comp_per_tissue_bar <- plot_grid(coloc_med_comp_per_tissue_bar, tgfm_med_comp_per_tissue_bar, ncol=1)
	output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_tgfm_coloc_mediated_component_per_tissue_bar_plot_", trait_name, ".pdf")
	ggsave(joint_med_comp_per_tissue_bar, file=output_file, width=7.2, height=6.0, units="in")


}



# Create Joint plot summarizing fraction of trait components that are non-mediated and mediated coming from each tissue
output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_susie_stacked_barplot_summary_joint.pdf")
joint_stacked_bar_plot_summary <- plot_grid(plotlist=stacked_barplot_joint_susie_list, ncol=2)
ggsave(joint_stacked_bar_plot_summary, file=output_file, width=9.2, height=9.4, units="in")


# Create joint plot summarizing fraction of trait components that are mediated and coming from each tissue
output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_susie_mediated_components_stacked_barplot_summary_joint.pdf")
joint_stacked_bar_plot_summary <- plot_grid(plotlist=stacked_barplot_joint_susie_mediated_list, ncol=2)
ggsave(joint_stacked_bar_plot_summary, file=output_file, width=9.2, height=9.4, units="in")


# Create plot summarizing fraction of trait components that are mediated in each trait
df <- data.frame(average_mediated_probability=avg_med_prob_arr, average_mediated_probability_se=se_med_prob_arr, trait=factor(trait_names_big,levels=trait_names_big))
output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_average_expression_mediated_probability_se_barplot.pdf")
med_prob_se_barplot <- make_mediated_prob_se_barplot(df)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=5.0, units="in")



}




# Load in gtex tissues
tissue_df <- read.table(gtex_pseudotissue_file,header=TRUE)
tissue_names <- as.character(tissue_df$pseudotissue_name)
tissue_names <- str_replace_all(tissue_names, "-", "_")

# Load in gtex tissue colors 
gtex_colors_df <- read.table(gtex_tissue_colors_file, header=TRUE, sep="\t")
gtex_colors_df$tissue_site_detail_id = as.character(gtex_colors_df$tissue_site_detail_id)
gtex_colors_df$tissue_site_detail_id[23] = "Cells_Cultured_fibroblasts"


# trait names
trait_df <- read.table(trait_file, header=TRUE)
trait_df = trait_df[as.character(trait_df$study_name) != "biochemistry_Triglycerides", ]
trait_df = trait_df[as.character(trait_df$study_name) != "mental_NEUROTICISM", ]
trait_df = trait_df[as.character(trait_df$study_name) != "disease_ALLERGY_ECZEMA_DIAGNOSED", ]

trait_names_big <- as.character(trait_df$study_name)
trait_names_big <- c("biochemistry_Cholesterol", "blood_WHITE_COUNT")
trait_names_big <- c("bp_DIASTOLICadjMEDz", "body_WHRadjBMIz")

# Model parameters
fusion_weights="False"
gene_version = "cis_heritable_genes"

stacked_barplot_joint_susie_mediated_list <- list()


# Loop through traits
for (itera in 1:length(trait_names_big)) {

	# Name of trait defining current loop
	trait_name <- trait_names_big[itera]

	# Get causal tissue and single_causal tissue corresponding to this trait
	causal_tissue_list <- get_causal_tissues_corresponding_to_trait(trait_name)
	causal_tissue <- causal_tissue_list[[1]]
	single_causal_tissue <- causal_tissue_list[[2]]

	# Load in organized TGFM results
	trait_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_component_organized_joint_tgfm_results.txt")
	trait_df <- read.table(trait_file, header=TRUE)
	

	# Calculate non mediated probability across all components for this trait
	susie_joint_non_med_prob_df <- extract_joint_susie_non_med_prob_df(trait_df, trait_name, tissue_names, "robust_tgfm_tissue_prior_joint_susie_prob")
	susie_joint_tissue_prior_df <- extract_joint_susie_df(trait_df, trait_name, tissue_names, "robust_tgfm_tissue_prior_joint_susie_prob")
	non_med_prob = 1.0 - (sum(susie_joint_tissue_prior_df$raw_statistic)/length(unique(susie_joint_tissue_prior_df$trait_component)))
	med_prob = 1.0 - non_med_prob

	# Summarize mediated and probabilities across components
	stacked_barplot_mediated_summary <- summarize_mediated_components_for_joint_susie_model_for_pres(susie_joint_tissue_prior_df, tissue_names, trait_name, gtex_colors_df) 
	stacked_barplot_joint_susie_mediated_list[[itera]] = stacked_barplot_mediated_summary

}



# Create joint plot summarizing fraction of trait components that are mediated and coming from each tissue
output_file <- paste0(rss_multivariate_twas_visualization_dir, "joint_susie_mediated_components_stacked_barplot_summary_joint_for_presentation.pdf")
joint_stacked_bar_plot_summary <- plot_grid(plotlist=stacked_barplot_joint_susie_mediated_list, ncol=1)
ggsave(joint_stacked_bar_plot_summary, file=output_file, width=14.5, height=5.5, units="in")



