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
trait_names <- as.character(trait_df$study_name)

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


