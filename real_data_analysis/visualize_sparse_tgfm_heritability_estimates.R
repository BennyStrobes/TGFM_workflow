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




get_number_of_genes_per_tissue <- function(num_genes_dir) {
	filer <- paste0(num_genes_dir, "number_of_reference_genes.txt")
	aa <- read.table(filer, header=TRUE,sep='\t')
	genes_per_tissue = aa$num_genes 


	return(genes_per_tissue)
}

get_numer_of_variants <- function(num_genes_and_variants_dir) {
	filer <- paste0(num_genes_and_variants_dir, "number_of_reference_variants.txt")
	aa = read.table(filer, header=FALSE)
	num_var = aa$V1[1]
	return(num_var)
}



load_in_gamma_ldsc_results <- function(tgfm_h2_results_dir, trait_names, num_predictors) {
	anno_name_arr <- c()
	trait_name_arr <- c()
	tau_arr <- c()
	tau_se_arr <- c()
	tau_z_arr <- c()
	h2_med_arr <- c()
	h2_med_se_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		# Load in summary file for trait
		trait_summary_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_top_window_False__remove_testis_True_results_organized.txt")
		trait_summary_df <- read.table(trait_summary_file, header=TRUE, sep="\t")

		# Remove intercept from trait summary df
		trait_summary_df_no_intercept = trait_summary_df[2:(dim(trait_summary_df)[1]),]

		anno_name_arr <- c(anno_name_arr, as.character(trait_summary_df_no_intercept$Annotation_name))
		trait_name_arr <- c(trait_name_arr, rep(trait_name, length(trait_summary_df_no_intercept$Annotation_name)))
		tau_arr <- c(tau_arr, trait_summary_df_no_intercept$tau)
		tau_se_arr <- c(tau_se_arr, trait_summary_df_no_intercept$tau_se)
		tau_z_arr <- c(tau_z_arr, trait_summary_df_no_intercept$tau_z)

		trait_h2_med = trait_summary_df_no_intercept$tau*num_predictors
		trait_h2_med_se = trait_summary_df_no_intercept$tau_se*num_predictors

		h2_med_arr <- c(h2_med_arr, trait_h2_med)
		h2_med_se_arr <- c(h2_med_se_arr, trait_h2_med_se)
	}


	# Put into neat data frame
	df <- data.frame(annotation=anno_name_arr, trait=trait_name_arr, tau=tau_arr, tau_se=tau_se_arr, tau_z=tau_z_arr, h2_med=h2_med_arr, h2_med_se=h2_med_se_arr)

	return(df)
}
load_in_sparse_gamma_ldsc_results <- function(sparse_h2_results_dir, trait_names, num_predictors) {
	anno_name_arr <- c()
	trait_name_arr <- c()
	tau_arr <- c()
	h2_med_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		# Load in summary file for trait
		trait_summary_file <- paste0(sparse_h2_results_dir, "sparse_ldsc_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_remove_testis_multivariate_results_sparse_loglss_sldsc_fixed_genotype_results_pmces.txt")
		trait_summary_df <- read.table(trait_summary_file, header=TRUE, sep="\t")


		anno_name_arr <- c(anno_name_arr, as.character(trait_summary_df$annotation))
		trait_name_arr <- c(trait_name_arr, rep(trait_name, length(trait_summary_df$annotation)))
		tau_arr <- c(tau_arr, trait_summary_df$pmces_tau)

		trait_h2_med = trait_summary_df$pmces_tau*num_predictors

		h2_med_arr <- c(h2_med_arr, trait_h2_med)
	}


	# Put into neat data frame
	df <- data.frame(annotation=anno_name_arr, trait=trait_name_arr, tau=tau_arr, tau_se=rep(0.0, length(tau_arr)), h2_med=h2_med_arr)

	return(df)
}


load_in_sparse_gaussian_ldsc_results <- function(sparse_h2_results_dir, trait_names, num_predictors) {
	anno_name_arr <- c()
	trait_name_arr <- c()
	tau_arr <- c()
	h2_med_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		# Load in summary file for trait
		trait_summary_file <- paste0(sparse_h2_results_dir, "sparse_ldsc_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_remove_testis_multivariate_results_sparse_sldsc_fixed_genotype_results_pmces.txt")
		trait_summary_df <- read.table(trait_summary_file, header=TRUE, sep="\t")


		anno_name_arr <- c(anno_name_arr, as.character(trait_summary_df$annotation))
		trait_name_arr <- c(trait_name_arr, rep(trait_name, length(trait_summary_df$annotation)))
		tau_arr <- c(tau_arr, trait_summary_df$pmces_tau)

		trait_h2_med = trait_summary_df$pmces_tau*num_predictors

		h2_med_arr <- c(h2_med_arr, trait_h2_med)
	}


	# Put into neat data frame
	df <- data.frame(annotation=anno_name_arr, trait=trait_name_arr, tau=tau_arr, tau_se=rep(0.0, length(tau_arr)), h2_med=h2_med_arr)

	return(df)
}


load_in_gaussian_ldsc_results <- function(sparse_h2_results_dir, trait_names, num_predictors) {
	anno_name_arr <- c()
	trait_name_arr <- c()
	tau_arr <- c()
	tau_se_arr <- c()
	tau_z_arr <- c()
	h2_med_arr <- c()
	h2_med_se_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		# Load in summary file for trait
		trait_summary_file <- paste0(sparse_h2_results_dir, "sparse_ldsc_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_remove_testis_multivariate_results_organized.txt")
		trait_summary_df <- read.table(trait_summary_file, header=TRUE, sep="\t")

		# Remove intercept from trait summary df
		trait_summary_df_no_intercept = trait_summary_df[2:(dim(trait_summary_df)[1]),]

		anno_name_arr <- c(anno_name_arr, as.character(trait_summary_df_no_intercept$Annotation_name))
		trait_name_arr <- c(trait_name_arr, rep(trait_name, length(trait_summary_df_no_intercept$Annotation_name)))
		tau_arr <- c(tau_arr, trait_summary_df_no_intercept$tau)
		tau_se_arr <- c(tau_se_arr, trait_summary_df_no_intercept$tau_se)
		tau_z_arr <- c(tau_z_arr, trait_summary_df_no_intercept$tau_z)

		trait_h2_med = trait_summary_df_no_intercept$tau*num_predictors
		trait_h2_med_se = trait_summary_df_no_intercept$tau_se*num_predictors

		h2_med_arr <- c(h2_med_arr, trait_h2_med)
		h2_med_se_arr <- c(h2_med_se_arr, trait_h2_med_se)
	}


	# Put into neat data frame
	df <- data.frame(annotation=anno_name_arr, trait=trait_name_arr, tau=tau_arr, tau_se=tau_se_arr, tau_z=tau_z_arr, h2_med=h2_med_arr, h2_med_se=h2_med_se_arr)

	return(df)
}
make_per_predictor_se_barplot2 <- function(tmp_gamma_ldsc_df, tmp_gaussian_ldsc_df, trait_name, model1_name, model2_name) {
	# organize data
	predictor_arr <- c(as.character(tmp_gamma_ldsc_df$annotation), as.character(tmp_gaussian_ldsc_df$annotation))
	tau_arr <- c(tmp_gamma_ldsc_df$tau, tmp_gaussian_ldsc_df$tau)
	tau_lb_arr <- c(tmp_gamma_ldsc_df$tau - 1.96*tmp_gamma_ldsc_df$tau_se, tmp_gaussian_ldsc_df$tau - 1.96*tmp_gaussian_ldsc_df$tau_se)
	tau_ub_arr <- c(tmp_gamma_ldsc_df$tau + 1.96*tmp_gamma_ldsc_df$tau_se, tmp_gaussian_ldsc_df$tau + 1.96*tmp_gaussian_ldsc_df$tau_se)
	method_arr <- c(rep(model1_name, length(tmp_gamma_ldsc_df$tau)), rep(model2_name, length(tmp_gaussian_ldsc_df$tau)))
	df <- data.frame(annotation=factor(predictor_arr), tau=tau_arr, tau_lb=tau_lb_arr, tau_ub=tau_ub_arr, method=method_arr)

	df$annotation = str_replace_all(as.character(df$annotation), "-", "_")
	df$annotation = str_replace_all(as.character(df$annotation), "Expression_", "")
	df$annotation <- recode(df$annotation, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	ordered_anno = as.character(df$annotation[df$method==model1_name])
	df$annotation = factor(df$annotation, levels=ordered_anno)



	p<-ggplot(data=df, aes(x=annotation, y=tau, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=tau_lb, ymax=tau_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=trait_name, x="") + theme(legend.position="top")

  	return(p)
}

make_per_predictor_se_barplot <- function(tmp_gamma_ldsc_df, tmp_gaussian_ldsc_df, trait_name) {
	# organize data
	predictor_arr <- c(as.character(tmp_gamma_ldsc_df$annotation), as.character(tmp_gaussian_ldsc_df$annotation))
	tau_arr <- c(tmp_gamma_ldsc_df$tau, tmp_gaussian_ldsc_df$tau)
	tau_lb_arr <- c(tmp_gamma_ldsc_df$tau - 1.96*tmp_gamma_ldsc_df$tau_se, tmp_gaussian_ldsc_df$tau - 1.96*tmp_gaussian_ldsc_df$tau_se)
	tau_ub_arr <- c(tmp_gamma_ldsc_df$tau + 1.96*tmp_gamma_ldsc_df$tau_se, tmp_gaussian_ldsc_df$tau + 1.96*tmp_gaussian_ldsc_df$tau_se)
	method_arr <- c(rep("SLDSC-Gamma", length(tmp_gamma_ldsc_df$tau)), rep("SLDSC-Gaussian", length(tmp_gaussian_ldsc_df$tau)))
	df <- data.frame(annotation=factor(predictor_arr), tau=tau_arr, tau_lb=tau_lb_arr, tau_ub=tau_ub_arr, method=method_arr)


	df$annotation = str_replace_all(as.character(df$annotation), "-", "_")
	df$annotation = str_replace_all(as.character(df$annotation), "Expression_", "")
	df$annotation <- recode(df$annotation, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	ordered_anno = as.character(df$annotation[df$method=="SLDSC-Gamma"])
	df$annotation = factor(df$annotation, levels=ordered_anno)

	p<-ggplot(data=df, aes(x=annotation, y=tau, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=tau_lb, ymax=tau_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=trait_name, x="")

  	return(p)
}

extract_gaussian_ldsc_fraction_of_h2_med_expression_df <- function(sparse_h2_results_dir, trait_names, num_ele_arr) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	fraction_med_arr <- c()
	fraction_med_se_arr <- c()
	total_h2_arr <- c()
	total_h2_arr_se <- c()

	n_var = num_ele_arr[1]
	n_genes_per_tiss = num_ele_arr[2:length(num_ele_arr)]


	for (trait_iter in 1:length(trait_names)) {

		trait_name <- trait_names[trait_iter]

		jacknifed_taus_file <- paste0(sparse_h2_results_dir, "sparse_ldsc_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_remove_testis_multivariate_results_raw_jacknifed_taus.txt")

		jacknifed_taus <- read.table(jacknifed_taus_file, header=FALSE, sep="\t")
		scaled_jacked_taus = jacknifed_taus*num_ele_arr

		genotype_h2 <- jacknifed_taus[,1]*n_var

		total_tiss_h2 <- rep(0, length(genotype_h2))

		for (tau_index in 2:(dim(scaled_jacked_taus)[2])) {
			total_tiss_h2 <- total_tiss_h2 + jacknifed_taus[,tau_index]*n_genes_per_tiss[(tau_index-1)]
		}


		fraction_h2_arr = total_tiss_h2/(total_tiss_h2+genotype_h2)
		fraction_mean = mean(fraction_h2_arr)
		diff_squared = (fraction_h2_arr - fraction_mean)**2
		num_jacknife_samples= length(fraction_h2_arr)
		jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
		jacknife_se = sqrt(jacknife_var)

		trait_names_arr <- c(trait_names_arr, trait_name)
		fraction_med_arr <- c(fraction_med_arr, fraction_mean)
		fraction_med_se_arr <- c(fraction_med_se_arr, jacknife_se)


		fraction_h2_arr = (total_tiss_h2+genotype_h2)
		fraction_mean = mean(fraction_h2_arr)
		diff_squared = (fraction_h2_arr - fraction_mean)**2
		num_jacknife_samples= length(fraction_h2_arr)
		jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
		jacknife_se = sqrt(jacknife_var)
		total_h2_arr <- c(total_h2_arr, fraction_mean)
		total_h2_arr_se <- c(total_h2_arr_se, jacknife_se)	
	}
	# Put all into compact df
	df <- data.frame(trait=trait_names_arr, fraction_h2_med=fraction_med_arr, fraction_h2_med_se=fraction_med_se_arr, total_h2=total_h2_arr, total_h2_se=total_h2_arr_se)
	return(df)
}


extract_ldsc_fraction_of_h2_med_expression_df <- function(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw, num_ele_arr) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	fraction_med_arr <- c()
	fraction_med_se_arr <- c()
	total_h2_arr <- c()
	total_h2_arr_se <- c()

	n_var = num_ele_arr[1]
	n_genes_per_tiss = num_ele_arr[2:length(num_ele_arr)]


	for (trait_iter in 1:length(trait_names)) {

		trait_name <- trait_names[trait_iter]

		# Load in jacknife data for this trait
		jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_top_window_False__remove_testis_True_jacknifed_mean_estimates.txt")
		jacknife_df <- read.table(jacknife_file, header=TRUE)

		tissue_name_raw <- "Genotype"
		genotype_df <- jacknife_df[as.character(jacknife_df$Class_name)==tissue_name_raw,]
		genotype_h2 <- genotype_df$h2*n_var

		total_tiss_h2 = rep(0.0, length(genotype_h2))



		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- tissue_names[tissue_iter]
			tissue_name_raw <- tissue_names_raw[tissue_iter]
			tissue_df <- jacknife_df[as.character(jacknife_df$Class_name)==tissue_name_raw,]
			tissue_h2 <- tissue_df$h2*n_genes_per_tiss[tissue_iter]
			total_tiss_h2 = total_tiss_h2 + tissue_h2
		}
		fraction_h2_arr = total_tiss_h2/(total_tiss_h2+genotype_h2)
		fraction_mean = mean(fraction_h2_arr)
		diff_squared = (fraction_h2_arr - fraction_mean)**2
		num_jacknife_samples= length(fraction_h2_arr)
		jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
		jacknife_se = sqrt(jacknife_var)

		trait_names_arr <- c(trait_names_arr, trait_name)
		fraction_med_arr <- c(fraction_med_arr, fraction_mean)
		fraction_med_se_arr <- c(fraction_med_se_arr, jacknife_se)


		fraction_h2_arr = (total_tiss_h2+genotype_h2)
		fraction_mean = mean(fraction_h2_arr)
		diff_squared = (fraction_h2_arr - fraction_mean)**2
		num_jacknife_samples= length(fraction_h2_arr)
		jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
		jacknife_se = sqrt(jacknife_var)
		total_h2_arr <- c(total_h2_arr, fraction_mean)
		total_h2_arr_se <- c(total_h2_arr_se, jacknife_se)	

	}

	# Put all into compact df
	df <- data.frame(trait=trait_names_arr, fraction_h2_med=fraction_med_arr, fraction_h2_med_se=fraction_med_se_arr, total_h2=total_h2_arr, total_h2_se=total_h2_arr_se)

	return(df)
}

make_fraction_mediated_h2_barplot <- function(gamma_ldsc_fraction_h2_df, gaussian_ldsc_fraction_h2_df) {
	frac_lb_arr <- c(gamma_ldsc_fraction_h2_df$fraction_h2_med - 1.96*gamma_ldsc_fraction_h2_df$fraction_h2_med_se, gaussian_ldsc_fraction_h2_df$fraction_h2_med - 1.96*gaussian_ldsc_fraction_h2_df$fraction_h2_med_se)
	frac_ub_arr <- c(gamma_ldsc_fraction_h2_df$fraction_h2_med + 1.96*gamma_ldsc_fraction_h2_df$fraction_h2_med_se, gaussian_ldsc_fraction_h2_df$fraction_h2_med + 1.96*gaussian_ldsc_fraction_h2_df$fraction_h2_med_se)
	method_arr <- c(rep("SLDSC-Gamma", length(gamma_ldsc_fraction_h2_df$fraction_h2_med)), rep("SLDSC-Gaussian", length(gaussian_ldsc_fraction_h2_df$fraction_h2_med)))
	frac_arr <- c(gamma_ldsc_fraction_h2_df$fraction_h2_med, gaussian_ldsc_fraction_h2_df$fraction_h2_med)
	trait_arr <- c(as.character(gamma_ldsc_fraction_h2_df$trait), as.character(gaussian_ldsc_fraction_h2_df$trait))
	df <- data.frame(trait=factor(trait_arr),fraction_h2=frac_arr, fraction_h2_lb=frac_lb_arr, fraction_h2_ub=frac_ub_arr, method=method_arr)

	p<-ggplot(data=df, aes(x=trait, y=fraction_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_h2_lb, ymax=fraction_h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title="", x="")

  	return(p)	

}
make_fraction_mediated_h2_barplot2 <- function(gamma_ldsc_fraction_h2_df, gaussian_ldsc_fraction_h2_df, model1_name, model2_name) {
	frac_lb_arr <- c(gamma_ldsc_fraction_h2_df$fraction_h2_med - 1.96*gamma_ldsc_fraction_h2_df$fraction_h2_med_se, gaussian_ldsc_fraction_h2_df$fraction_h2_med - 1.96*gaussian_ldsc_fraction_h2_df$fraction_h2_med_se)
	frac_ub_arr <- c(gamma_ldsc_fraction_h2_df$fraction_h2_med + 1.96*gamma_ldsc_fraction_h2_df$fraction_h2_med_se, gaussian_ldsc_fraction_h2_df$fraction_h2_med + 1.96*gaussian_ldsc_fraction_h2_df$fraction_h2_med_se)
	method_arr <- c(rep(model1_name, length(gamma_ldsc_fraction_h2_df$fraction_h2_med)), rep(model2_name, length(gaussian_ldsc_fraction_h2_df$fraction_h2_med)))
	frac_arr <- c(gamma_ldsc_fraction_h2_df$fraction_h2_med, gaussian_ldsc_fraction_h2_df$fraction_h2_med)
	trait_arr <- c(as.character(gamma_ldsc_fraction_h2_df$trait), as.character(gaussian_ldsc_fraction_h2_df$trait))
	df <- data.frame(trait=factor(trait_arr),fraction_h2=frac_arr, fraction_h2_lb=frac_lb_arr, fraction_h2_ub=frac_ub_arr, method=method_arr)

	p<-ggplot(data=df, aes(x=trait, y=fraction_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_h2_lb, ymax=fraction_h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title="", x="")

  	return(p)	

}

make_h2_barplot2 <- function(gamma_ldsc_fraction_h2_df, gaussian_ldsc_fraction_h2_df, model1_name, model2_name) {
	frac_lb_arr <- c(gamma_ldsc_fraction_h2_df$total_h2 - 1.96*gamma_ldsc_fraction_h2_df$total_h2_se, gaussian_ldsc_fraction_h2_df$total_h2 - 1.96*gaussian_ldsc_fraction_h2_df$total_h2_se)
	frac_ub_arr <- c(gamma_ldsc_fraction_h2_df$total_h2 + 1.96*gamma_ldsc_fraction_h2_df$total_h2_se, gaussian_ldsc_fraction_h2_df$total_h2 + 1.96*gaussian_ldsc_fraction_h2_df$total_h2_se)
	method_arr <- c(rep(model1_name, length(gamma_ldsc_fraction_h2_df$fraction_h2_med)), rep(model2_name, length(gaussian_ldsc_fraction_h2_df$fraction_h2_med)))
	frac_arr <- c(gamma_ldsc_fraction_h2_df$total_h2, gaussian_ldsc_fraction_h2_df$total_h2)
	trait_arr <- c(as.character(gamma_ldsc_fraction_h2_df$trait), as.character(gaussian_ldsc_fraction_h2_df$trait))
	df <- data.frame(trait=factor(trait_arr),fraction_h2=frac_arr, fraction_h2_lb=frac_lb_arr, fraction_h2_ub=frac_ub_arr, method=method_arr)

	p<-ggplot(data=df, aes(x=trait, y=fraction_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_h2_lb, ymax=fraction_h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title="", x="", y="total h2")

  	return(p)	

}

make_h2_barplot <- function(gamma_ldsc_fraction_h2_df, gaussian_ldsc_fraction_h2_df) {
	frac_lb_arr <- c(gamma_ldsc_fraction_h2_df$total_h2 - 1.96*gamma_ldsc_fraction_h2_df$total_h2_se, gaussian_ldsc_fraction_h2_df$total_h2 - 1.96*gaussian_ldsc_fraction_h2_df$total_h2_se)
	frac_ub_arr <- c(gamma_ldsc_fraction_h2_df$total_h2 + 1.96*gamma_ldsc_fraction_h2_df$total_h2_se, gaussian_ldsc_fraction_h2_df$total_h2 + 1.96*gaussian_ldsc_fraction_h2_df$total_h2_se)
	method_arr <- c(rep("SLDSC-Gamma", length(gamma_ldsc_fraction_h2_df$fraction_h2_med)), rep("SLDSC-Gaussian", length(gaussian_ldsc_fraction_h2_df$fraction_h2_med)))
	frac_arr <- c(gamma_ldsc_fraction_h2_df$total_h2, gaussian_ldsc_fraction_h2_df$total_h2)
	trait_arr <- c(as.character(gamma_ldsc_fraction_h2_df$trait), as.character(gaussian_ldsc_fraction_h2_df$trait))
	df <- data.frame(trait=factor(trait_arr),fraction_h2=frac_arr, fraction_h2_lb=frac_lb_arr, fraction_h2_ub=frac_ub_arr, method=method_arr)

	p<-ggplot(data=df, aes(x=trait, y=fraction_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_h2_lb, ymax=fraction_h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title="", x="", y="total h2")

  	return(p)	

}

make_tau_heatmap <- function(gamma_ldsc_df) {
	df <- gamma_ldsc_df[gamma_ldsc_df$annotation != "Genotype",]


	df$annotation = str_replace_all(as.character(df$annotation), "Expression_", "")
	df$annotation = str_replace_all(as.character(df$annotation), "-", "_")
	df$annotation <- recode(df$annotation, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	df$trait <- recode(df$trait, blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="blood_Reticulocyte_count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="blood_Hemoglobin")
	#df$tau_z[df$tau_z < threshold] = 0.0
	df$tau[df$tau < 0.0] = 0.0
	p <- ggplot(df, aes(x = annotation, y = trait, fill = tau)) +
  		geom_tile() +
  		theme(text = element_text(size=11), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="Expression-mediated\ntrait heritability tau",x="Tissue", y="GWAS trait",title="")

  	return(p)


}

make_z_score_heatmap <- function(gamma_ldsc_df, threshold=2.0) {
	df <- gamma_ldsc_df[gamma_ldsc_df$annotation != "Genotype",]


	df$annotation = str_replace_all(as.character(df$annotation), "Expression_", "")
	df$annotation = str_replace_all(as.character(df$annotation), "-", "_")
	df$annotation <- recode(df$annotation, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	df$trait <- recode(df$trait, blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="blood_Reticulocyte_count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="blood_Hemoglobin")
	#df$tau_z[df$tau_z < threshold] = 0.0
	df$tau_z = abs(df$tau_z)
	p <- ggplot(df, aes(x = annotation, y = trait, fill = tau_z)) +
  		geom_tile() +
  		theme(text = element_text(size=11), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="Expression-mediated\ntrait heritability z-score",x="Tissue", y="GWAS trait",title="")

  	return(p)


}


get_fraction_h2_for_sparse_models <- function(trait_names, sparse_gaussian_ldsc_df, num_) {

	trait_names_arr <- c()
	fraction_med_arr <- c()
	total_h2_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		tmp_df <- sparse_gaussian_ldsc_df[sparse_gaussian_ldsc_df$trait==trait_name,]
		geno_h2 = sum(tmp_df$h2_med[tmp_df$annotation=="Genotype"])
		expr_h2 = sum(tmp_df$h2_med[tmp_df$annotation!="Genotype"])

		trait_names_arr <- c(trait_names_arr, trait_name)
		fraction_med_arr <- c(fraction_med_arr, expr_h2/(expr_h2+geno_h2))
		total_h2_arr <- c(total_h2_arr, (expr_h2+geno_h2))

	}


	df <- data.frame(trait=trait_names_arr, fraction_h2_med=fraction_med_arr, fraction_h2_med_se=rep(0.0, length(fraction_med_arr)), total_h2=total_h2_arr, total_h2_se=rep(0.0, length(fraction_med_arr)))


}

######################
# Command line args
#######################
tissue_names_file = args[1]
trait_names_file = args[2]
num_genes_and_variants_dir = args[3]
tgfm_h2_results_dir = args[4]
sparse_h2_results_dir = args[5]
viz_dir = args[6]


# Load in gtex tissues
tissue_df <- read.table(tissue_names_file,header=TRUE)
tissue_names_raw <- as.character(tissue_df$pseudotissue_name)
tissue_names <- str_replace_all(tissue_names_raw, "-", "_")
valid_tiss = tissue_names != "Testis"
num_genes_per_tissue <- get_number_of_genes_per_tissue(num_genes_and_variants_dir)

num_genes_per_tissue = num_genes_per_tissue[valid_tiss]
tissue_names = tissue_names[valid_tiss]
tissue_names_raw = tissue_names_raw[valid_tiss]

num_variants <- get_numer_of_variants(num_genes_and_variants_dir)
num_features = c(num_variants, num_genes_per_tissue)


# Get list of trait names
trait_names_df <- read.table(trait_names_file, header=TRUE,sep="\t")
trait_names <- as.character(trait_names_df$study_name)




####################################################################
# First compare gamma-sldsc regression vs gaussian-sldsc regression
####################################################################
# Load in per feature data
gamma_ldsc_df <- load_in_gamma_ldsc_results(tgfm_h2_results_dir, trait_names, num_features)
gaussian_ldsc_df <- load_in_gaussian_ldsc_results(sparse_h2_results_dir, trait_names, num_features)
# Load in aggregated data
gamma_ldsc_fraction_h2_df <- extract_ldsc_fraction_of_h2_med_expression_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw, num_features)
gaussian_ldsc_fraction_h2_df <- extract_gaussian_ldsc_fraction_of_h2_med_expression_df(sparse_h2_results_dir, trait_names, num_features)

if (FALSE) {
# Make fraction mediated h2 barplot for two models
fraction_mediated_barplot <- make_fraction_mediated_h2_barplot(gamma_ldsc_fraction_h2_df, gaussian_ldsc_fraction_h2_df)
output_file <- paste0(viz_dir, "fraction_h2_mediated_se_barplot.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")

# Make h2 barplot for two models
fraction_mediated_barplot <- make_h2_barplot(gamma_ldsc_fraction_h2_df, gaussian_ldsc_fraction_h2_df)
output_file <- paste0(viz_dir, "h2_se_barplot.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")


# Make z-score boxplot based on both
z_score_heatmap <- make_z_score_heatmap(gamma_ldsc_df, threshold=0.0)
output_file <- paste0(viz_dir, "gamma_ldsc_tau_z_score_heatmap.pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=6.0, units="in")

z_score_heatmap <- make_z_score_heatmap(gaussian_ldsc_df, threshold=0.0)
output_file <- paste0(viz_dir, "gaussian_ldsc_tau_z_score_heatmap.pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=6.0, units="in")


for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]

	# Make se bar plot show per predictor h2
	tmp_gamma_ldsc_df = gamma_ldsc_df[gamma_ldsc_df$trait==trait_name, ]
	tmp_gaussian_ldsc_df = gaussian_ldsc_df[gaussian_ldsc_df$trait==trait_name, ]
	predictor_se_barplot <- make_per_predictor_se_barplot(tmp_gamma_ldsc_df, tmp_gaussian_ldsc_df, trait_name)
	output_file <- paste0(viz_dir, "trait_specific_", trait_name, "_per_predictor_tau_se_barplot.pdf")
	ggsave(predictor_se_barplot, file=output_file, width=7.2, height=4.0, units="in")
}

}




###########################################
# Now compare to sparse regression
###########################################
sparse_gaussian_ldsc_df <- load_in_sparse_gaussian_ldsc_results(sparse_h2_results_dir, trait_names, num_features)
sparse_gaussian_ldsc_fraction_h2_df <- get_fraction_h2_for_sparse_models(trait_names, sparse_gaussian_ldsc_df, num_features)
sparse_gamma_ldsc_df <- load_in_sparse_gamma_ldsc_results(sparse_h2_results_dir, trait_names, num_features)
sparse_gamma_ldsc_fraction_h2_df <- get_fraction_h2_for_sparse_models(trait_names, sparse_gamma_ldsc_df, num_features)




sparse_tau_heatmap <- make_tau_heatmap(sparse_gaussian_ldsc_df)
output_file <- paste0(viz_dir, "gaussian_sparse_ldsc_tau_heatmap.pdf")
ggsave(sparse_tau_heatmap, file=output_file, width=7.2, height=6.0, units="in")


sparse_tau_heatmap <- make_tau_heatmap(sparse_gamma_ldsc_df)
output_file <- paste0(viz_dir, "gamma_sparse_ldsc_tau_heatmap.pdf")
ggsave(sparse_tau_heatmap, file=output_file, width=7.2, height=6.0, units="in")

# Make fraction mediated h2 barplot for two models
fraction_mediated_barplot <- make_fraction_mediated_h2_barplot2(gaussian_ldsc_fraction_h2_df, sparse_gaussian_ldsc_fraction_h2_df, "SLDSC-Gaussian", "Sparse-SLDSC-Gaussian")
output_file <- paste0(viz_dir, "fraction_h2_mediated_se_barplot_sparse_model.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")

fraction_mediated_barplot <- make_fraction_mediated_h2_barplot2(gamma_ldsc_fraction_h2_df, sparse_gamma_ldsc_fraction_h2_df, "SLDSC-Gamma", "Sparse-SLDSC-Gamma")
output_file <- paste0(viz_dir, "fraction_h2_mediated_se_barplot_gamma_sparse_model.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")


# Make h2 barplot for two models
fraction_mediated_barplot <- make_h2_barplot2(gaussian_ldsc_fraction_h2_df, sparse_gaussian_ldsc_fraction_h2_df, "SLDSC-Gaussian", "Sparse-SLDSC-Gaussian")
output_file <- paste0(viz_dir, "h2_se_barplot_sparse_model.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")

fraction_mediated_barplot <- make_h2_barplot2(gamma_ldsc_fraction_h2_df, sparse_gamma_ldsc_fraction_h2_df, "SLDSC-Gamma", "Sparse-SLDSC-Gamma")
output_file <- paste0(viz_dir, "h2_se_barplot_gamma_sparse_model.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")


for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	# Make se bar plot show per predictor h2
	tmp_gaussian_ldsc_df = gaussian_ldsc_df[gaussian_ldsc_df$trait==trait_name, ]
	tmp_sparse_gaussian_ldsc_df = sparse_gaussian_ldsc_df[sparse_gaussian_ldsc_df$trait==trait_name, ]
	predictor_se_barplot <- make_per_predictor_se_barplot2(tmp_gaussian_ldsc_df, tmp_sparse_gaussian_ldsc_df, trait_name, "SLDSC-Gaussian", "Sparse-SLDSC-Gaussian")
	output_file <- paste0(viz_dir, "trait_specific_", trait_name, "_per_predictor_sparse_tau_se_barplot.pdf")
	ggsave(predictor_se_barplot, file=output_file, width=7.2, height=4.0, units="in")

	# Make se bar plot show per predictor h2
	tmp_gamma_ldsc_df = gamma_ldsc_df[gamma_ldsc_df$trait==trait_name, ]
	tmp_sparse_gamma_ldsc_df = sparse_gamma_ldsc_df[sparse_gamma_ldsc_df$trait==trait_name, ]
	predictor_se_barplot <- make_per_predictor_se_barplot2(tmp_gamma_ldsc_df, tmp_sparse_gamma_ldsc_df, trait_name, "SLDSC-Gamma", "Sparse-SLDSC-Gamma")
	output_file <- paste0(viz_dir, "trait_specific_", trait_name, "_per_predictor_gamma_sparse_tau_se_barplot.pdf")
	ggsave(predictor_se_barplot, file=output_file, width=7.2, height=4.0, units="in")

}

