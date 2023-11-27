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

make_rss_variance_proportion_heatmap <- function(df, ordered_trait_names) {
	variance_proportion_arr <- c()
	tissue_names_arr <- c()
	trait_names_arr <- c()

	print(head(df))

	for (trait_iter in 1:length(ordered_trait_names)) {
		trait_name <- ordered_trait_names[trait_iter]
		df_sub = df[as.character(df$trait) == trait_name,]
		sub_tissue_arr = as.character(df_sub$tissue)
		sub_variance_prop = df_sub$per_gene_h2/sum(df_sub$per_gene_h2)

		variance_proportion_arr <- c(variance_proportion_arr, sub_variance_prop)
		tissue_names_arr <- c(tissue_names_arr, sub_tissue_arr)
		trait_names_arr <- c(trait_names_arr, rep(trait_name, length(sub_variance_prop)))
	}



	# PUt in inorganized data frame
	df2 <- data.frame(trait=as.character(trait_names_arr), tissue=as.character(tissue_names_arr), variance_proportion=variance_proportion_arr)
	df2$trait = factor(df2$trait, levels=ordered_trait_names)

	#df2$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	#df2$tissue <- recode(df2$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	#df2$trait <- recode(df2$trait, biochemistry_Cholesterol="Cholesterol", blood_EOSINOPHIL_COUNT="Blood eosinophil count", blood_RBC_DISTRIB_WIDTH="Blood RBC width", blood_RED_COUNT="Blood red count", blood_WHITE_COUNT="Blood white count", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Corposcular hemoglobin", disease_HYPERTENSION_DIAGNOSED="Hypertension")
	#df$trait <- factor(df$trait, levels=c("Eczema", "College Education", "Diastolic BP", "WHR-adj BMI", "Blood white count", "Blood red count", "Blood RBC width", "Blood eosinophil count", "Cholesterol"))
	p <- ggplot(df2, aes(x = tissue, y = trait, fill = variance_proportion)) +
  		geom_tile() +
  		theme(text = element_text(size=10), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="Proportion of expression-mediated\ntrait heritability",x="Tissue", y="GWAS trait",title="")
  	return(p)
}


make_variance_proportion_heatmap <- function(tgfm_h2_results_dir, trait_names, tissue_names) {
	variance_proportion_arr <- c()
	tissue_names_arr <- c()
	trait_names_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		# File name for trait
		variance_file = paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_learn_intercept_mean_estimates.txt")
		# load in df
		variance_df <- read.table(variance_file, header=TRUE)
		# filter to genes
		variance_df = variance_df[3:(dim(variance_df)[1]),]

		trait_tissue_names <- as.character(variance_df$Class_name)
		trait_var = variance_df$h2
		if (max(trait_var) > 1e-6) {
		norm_trait_var = trait_var/sum(trait_var)
		
		# Add to normalized arr
		variance_proportion_arr <- c(variance_proportion_arr, norm_trait_var)
		tissue_names_arr <- c(tissue_names_arr, trait_tissue_names)
		trait_names_arr <- c(trait_names_arr, rep(trait_name, length(trait_tissue_names)))
		}
	}

	# PUt in inorganized data frame
	df <- data.frame(variance_proportion=variance_proportion_arr, tissue=tissue_names_arr, trait=trait_names_arr)
	df$trait = factor(df$trait, levels=trait_names)

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	#df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", blood_EOSINOPHIL_COUNT="Blood eosinophil count", blood_RBC_DISTRIB_WIDTH="Blood RBC width", blood_RED_COUNT="Blood red count", blood_WHITE_COUNT="Blood white count", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema")
	#df$trait <- factor(df$trait, levels=c("Eczema", "College Education", "Diastolic BP", "WHR-adj BMI", "Blood white count", "Blood red count", "Blood RBC width", "Blood eosinophil count", "Cholesterol"))
	p <- ggplot(df, aes(x = tissue, y = trait, fill = variance_proportion)) +
  		geom_tile() +
  		theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="Proportion of expression-mediated\ntrait heritability",x="Tissue", y="GWAS trait",title="")
  	return(p)
}

make_z_score_heatmap <- function(tgfm_h2_results_dir, trait_names, tissue_names) {
	variance_proportion_arr <- c()
	tissue_names_arr <- c()
	trait_names_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		# File name for trait
		variance_file = paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_top_window_False_mean_estimates.txt")
		# load in df
		variance_df <- read.table(variance_file, header=TRUE)
		# filter to genes
		variance_df = variance_df[3:(dim(variance_df)[1]),]

		trait_tissue_names <- as.character(variance_df$Class_name)
		trait_var = variance_df$h2
		
		jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_top_window_False_jacknifed_mean_estimates.txt")
		jacknife_df <- read.table(jacknife_file, header=TRUE)
		z_score_vec = rep(0.0, length(trait_var))
		if (max(trait_var) > 1e-6) {
		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- tissue_names[tissue_iter]
			tissue_df <- jacknife_df[as.character(jacknife_df$Class_name)==tissue_name,]
			num_jacknife_samples = dim(tissue_df)[1]

			h2_mean = mean(tissue_df$h2)
			diff_squared = (tissue_df$h2 - h2_mean)**2
			jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
			jacknife_se = sqrt(jacknife_var)
			h2_lb = h2_mean - 1.96*jacknife_se
			z_score = h2_mean/jacknife_se
			if (h2_lb < 0.0) {
				z_score = 0.0
			}
			
			z_score_vec[tissue_iter] = z_score
		}



		# Add to normalized arr
		variance_proportion_arr <- c(variance_proportion_arr, z_score_vec)
		tissue_names_arr <- c(tissue_names_arr, trait_tissue_names)
		trait_names_arr <- c(trait_names_arr, rep(trait_name, length(trait_tissue_names)))
		}
	}

	# PUt in inorganized data frame
	df <- data.frame(variance_proportion=variance_proportion_arr, tissue=tissue_names_arr, trait=trait_names_arr)
	df$trait = factor(df$trait, levels=trait_names)

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	#df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", blood_EOSINOPHIL_COUNT="Blood eosinophil count", blood_RBC_DISTRIB_WIDTH="Blood RBC width", blood_RED_COUNT="Blood red count", blood_WHITE_COUNT="Blood white count", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema")
	#df$trait <- factor(df$trait, levels=c("Eczema", "College Education", "Diastolic BP", "WHR-adj BMI", "Blood white count", "Blood red count", "Blood RBC width", "Blood eosinophil count", "Cholesterol"))
	p <- ggplot(df, aes(x = tissue, y = trait, fill = variance_proportion)) +
  		geom_tile() +
  		theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="Expression-mediated\ntrait heritability z-score",x="Tissue", y="GWAS trait",title="")
  	return(p)
}


jacknife_mean_standard_error_plot_for_single_trait <- function(trait_name, tgfm_h2_results_dir, tissue_names) {
	jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_learn_intercept_jacknifed_mean_estimates.txt")
	jacknife_df <- read.table(jacknife_file, header=TRUE)

	# Initialize arrays to save data
	tissue_names_arr <- c()
	per_gene_h2_arr <- c()
	per_gene_h2_ub_arr <- c()
	per_gene_h2_lb_arr <- c()
	significant_arr <- c()
	color_arr <- c()


	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		tissue_df <- jacknife_df[as.character(jacknife_df$Class_name)==tissue_name,]
		num_jacknife_samples = dim(tissue_df)[1]

		h2_mean = mean(tissue_df$h2)
		diff_squared = (tissue_df$h2 - h2_mean)**2
		jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
		jacknife_se = sqrt(jacknife_var)

		h2_ub = h2_mean + 1.96*jacknife_se 
		h2_lb = h2_mean - 1.96*jacknife_se
		#print(tissue_name)
		#print(paste0("[", h2_lb, " , ", h2_ub, "]"))
		z_score = h2_mean/jacknife_se
		#print(z_score)
		if (h2_lb <= 0) {
			significant_arr <- c(significant_arr, "not_significant")
			color_arr <- c(color_arr, "#56B4E9")
		} else {
			significant_arr <- c(significant_arr, "significant")
			color_arr <- c(color_arr, "#099999")
		}

		tissue_names_arr <- c(tissue_names_arr, tissue_name)
		per_gene_h2_arr <- c(per_gene_h2_arr, h2_mean)
		per_gene_h2_lb_arr <- c(per_gene_h2_lb_arr, h2_lb)
		per_gene_h2_ub_arr <- c(per_gene_h2_ub_arr, h2_ub)

	}

	df <- data.frame(tissue=tissue_names_arr, gene_h2=per_gene_h2_arr, gene_h2_lb=per_gene_h2_lb_arr, gene_h2_ub=per_gene_h2_ub_arr, sig=significant_arr)

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	p<- ggplot(df, aes(x=tissue, y=gene_h2, fill=sig)) + 
 		 geom_bar(stat="identity", color="black", position=position_dodge()) +
  		 scale_fill_manual(values=c("#999999", "#56B4E9")) +
 		 figure_theme() +
  		 geom_errorbar(aes(ymin=gene_h2_lb, ymax=gene_h2_ub), width=.2, position=position_dodge(.9)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		 theme(legend.position="none")+
  		 labs(x="", y="Per gene h2", title=trait_name)
  	return(p)
}

extract_ldsc_per_gene_h2_plus_non_mediated_df <- function(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	tissue_names_arr <- c()
	per_gene_h2_arr <- c()
	per_gene_h2_se_arr <- c()	


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		# Load in jacknife data for this trait
		jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_top_window_False_jacknifed_mean_estimates.txt")
		jacknife_df <- read.table(jacknife_file, header=TRUE)

		tissue_name_raw <- "Genotype"
		tissue_name <- "Non-mediated"
		tissue_df <- jacknife_df[as.character(jacknife_df$Class_name)==tissue_name_raw,]
		num_jacknife_samples = dim(tissue_df)[1]

		h2_mean = mean(tissue_df$h2)
		diff_squared = (tissue_df$h2 - h2_mean)**2
		jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
		jacknife_se = sqrt(jacknife_var)


		trait_names_arr <- c(trait_names_arr, trait_name)
		tissue_names_arr <- c(tissue_names_arr, tissue_name)
		per_gene_h2_arr <- c(per_gene_h2_arr, h2_mean)
		per_gene_h2_se_arr <- c(per_gene_h2_se_arr, jacknife_se)	


		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- tissue_names[tissue_iter]
			tissue_name_raw <- tissue_names_raw[tissue_iter]
			tissue_df <- jacknife_df[as.character(jacknife_df$Class_name)==tissue_name_raw,]
			num_jacknife_samples = dim(tissue_df)[1]

			h2_mean = mean(tissue_df$h2)
			diff_squared = (tissue_df$h2 - h2_mean)**2
			jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
			jacknife_se = sqrt(jacknife_var)


			trait_names_arr <- c(trait_names_arr, trait_name)
			tissue_names_arr <- c(tissue_names_arr, tissue_name)
			per_gene_h2_arr <- c(per_gene_h2_arr, h2_mean)
			per_gene_h2_se_arr <- c(per_gene_h2_se_arr, jacknife_se)	
		}
	}

	# Put all into compact df
	df <- data.frame(trait=trait_names_arr, tissue=tissue_names_arr, per_gene_h2=per_gene_h2_arr, per_gene_h2_se=per_gene_h2_se_arr)

	return(df)

}

extract_ldsc_per_gene_h2_df <- function(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	tissue_names_arr <- c()
	per_gene_h2_arr <- c()
	per_gene_h2_se_arr <- c()	


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		# Load in jacknife data for this trait
		jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_top_window_False_jacknifed_mean_estimates.txt")
		jacknife_df <- read.table(jacknife_file, header=TRUE)

		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- tissue_names[tissue_iter]
			tissue_name_raw <- tissue_names_raw[tissue_iter]
			tissue_df <- jacknife_df[as.character(jacknife_df$Class_name)==tissue_name_raw,]
			num_jacknife_samples = dim(tissue_df)[1]

			h2_mean = mean(tissue_df$h2)
			diff_squared = (tissue_df$h2 - h2_mean)**2
			jacknife_var = sum(diff_squared)*(num_jacknife_samples-1.0)/(num_jacknife_samples)
			jacknife_se = sqrt(jacknife_var)


			trait_names_arr <- c(trait_names_arr, trait_name)
			tissue_names_arr <- c(tissue_names_arr, tissue_name)
			per_gene_h2_arr <- c(per_gene_h2_arr, h2_mean)
			per_gene_h2_se_arr <- c(per_gene_h2_se_arr, jacknife_se)	
		}
	}

	# Put all into compact df
	df <- data.frame(trait=trait_names_arr, tissue=tissue_names_arr, per_gene_h2=per_gene_h2_arr, per_gene_h2_se=per_gene_h2_se_arr)

	return(df)

}

extract_ldsc_fraction_of_h2_med_expression_df2<- function(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw, num_ele_arr) {


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
		jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_top_window_True_jacknifed_mean_estimates.txt")
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
		jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_cis_heritable_gene_learn_intercept_jacknifed_mean_estimates.txt")
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


extract_rss_fraction_of_h2_med_expression_df <- function(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw, num_ele_arr) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	fraction_med_arr <- c()
	total_h2_arr <- c()

	n_var = num_ele_arr[1]
	n_genes_per_tiss = num_ele_arr[2:length(num_ele_arr)]


	for (trait_iter in 1:length(trait_names)) {

		trait_name <- trait_names[trait_iter]

		non_med_rss_file <- paste0(tgfm_h2_results_dir, "tgfm_rss_likelihood_parallel_style_heritability_", trait_name, "_cis_heritable_gene_standardize_expr_True__robust_pleiotropic_prior_precision_temp.txt")
		non_med_rss_df <- read.table(non_med_rss_file, header=TRUE)

		geno_h2_mean = non_med_rss_df$gamma_beta_b/non_med_rss_df$gamma_beta_a
		geno_h2 = geno_h2_mean*n_var

		rss_file <- paste0(tgfm_h2_results_dir, "tgfm_rss_likelihood_parallel_style_heritability_", trait_name, "_cis_heritable_gene_standardize_expr_True_robust_tissue_specific_prior_precision_temp.txt")
		rss_df <- read.table(rss_file, header=TRUE)

		gene_h2 <- rss_df$gamma_alpha_b/rss_df$gamma_alpha_a
		#print(gene_h2)
		total_gene_h2 <- sum(gene_h2*n_genes_per_tiss)

		trait_names_arr <- c(trait_names_arr, trait_name)
		fraction_med_arr <- c(fraction_med_arr, total_gene_h2/(total_gene_h2+geno_h2))

		total_h2_arr <- c(total_h2_arr, (total_gene_h2+geno_h2))

	}

	df <- data.frame(trait=trait_names_arr, fraction_h2_med=fraction_med_arr, total_h2=total_h2_arr)
	return(df)
}

extract_standard_sldsc_non_eqtl_cpp_removed_fraction_of_h2_med_expression_df <- function(standard_sldsc_results_dir, sldsc_anno_dir, trait_names, tissue_names, tissue_names_raw, suffix) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	fraction_med_arr <- c()
	total_h2_arr <- c()


	for (chrom_num in 1:22) {
		file_name <- paste0(sldsc_anno_dir, "baselineld_non_eqtl_cpp_removed.", chrom_num, ".", suffix)
		data = read.table(file_name, header=FALSE)
		if (chrom_num==1) {
			nvar = data[1,]
		} else {
			nvar = nvar + data[1,]
		}
	}

	for (trait_iter in 1:length(trait_names)) {

		trait_name <- trait_names[trait_iter]

		trait_file <- paste0(standard_sldsc_results_dir, "UKB_460K.", trait_name, "_sldsc_res_baselineld_non_eqtl_cpp_removed_.results")
		trait_df <- read.table(trait_file, header=TRUE, sep="\t")
		tau = trait_df$Coefficient
		per_anno_h2 = tau*nvar
		expr = per_anno_h2[73]
		expr = as.numeric(expr[1])[1]

		tot = sum(per_anno_h2)[1]



		total_h2_arr <- c(total_h2_arr, tot)
		fraction_med_arr <- c(fraction_med_arr, expr/(tot))
		trait_names_arr <- c(trait_names_arr, trait_name)
	}	
	df <- data.frame(trait=trait_names_arr, fraction_h2_med=fraction_med_arr, total_h2=total_h2_arr)
	return(df)
}


extract_standard_sldsc_all_anno_fraction_of_h2_med_expression_df <- function(standard_sldsc_results_dir, sldsc_anno_dir, trait_names, tissue_names, tissue_names_raw, suffix) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	fraction_med_arr <- c()
	total_h2_arr <- c()


	for (chrom_num in 1:22) {
		file_name <- paste0(sldsc_anno_dir, "baselineLD.", chrom_num, ".", suffix)
		data = read.table(file_name, header=FALSE)
		if (chrom_num==1) {
			nvar = data[1,]
		} else {
			nvar = nvar + data[1,]
		}
	}

	for (trait_iter in 1:length(trait_names)) {

		trait_name <- trait_names[trait_iter]

		trait_file <- paste0(standard_sldsc_results_dir, "UKB_460K.", trait_name, "_sldsc_res_.results")
		trait_df <- read.table(trait_file, header=TRUE, sep="\t")
		tau = trait_df$Coefficient
		per_anno_h2 = tau*nvar
		expr = per_anno_h2[73]
		expr = as.numeric(expr[1])[1]

		tot = sum(per_anno_h2)[1]



		total_h2_arr <- c(total_h2_arr, tot)
		fraction_med_arr <- c(fraction_med_arr, expr/(tot))
		trait_names_arr <- c(trait_names_arr, trait_name)
	}	
	df <- data.frame(trait=trait_names_arr, fraction_h2_med=fraction_med_arr, total_h2=total_h2_arr)
	return(df)
}


extract_standard_sldsc_fraction_of_h2_med_expression_df <- function(standard_sldsc_results_dir, standard_sldsc_processed_data_dir, trait_names, tissue_names, tissue_names_raw, suffix) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	fraction_med_arr <- c()
	total_h2_arr <- c()

	nvar1 = 0
	nvar2 = 0
	for (chrom_num in 1:22) {
		file_name <- paste0(standard_sldsc_processed_data_dir, "baseline_intercept_only.", chrom_num, ".", suffix)
		data = read.table(file_name, header=FALSE)
		nvar1 = nvar1 + data$V1[1]
		nvar2 = nvar2 + data$V2[1]
	}

	for (trait_iter in 1:length(trait_names)) {

		trait_name <- trait_names[trait_iter]

		trait_file <- paste0(standard_sldsc_results_dir, "UKB_460K.", trait_name, "_sldsc_res_intercept_only_.results")
		trait_df <- read.table(trait_file, header=TRUE, sep="\t")
		tau = trait_df$Coefficient
		tau_se = trait_df$Coefficient_std_error

		geno = tau[1]*nvar1
		expr = tau[2]*nvar2

		total_h2_arr <- c(total_h2_arr, geno+expr)
		fraction_med_arr <- c(fraction_med_arr, expr/(geno+expr))
		trait_names_arr <- c(trait_names_arr, trait_name)
	}	
	df <- data.frame(trait=trait_names_arr, fraction_h2_med=fraction_med_arr, total_h2=total_h2_arr)
	return(df)
}

extract_rss_per_gene_h2_plus_non_med_df <- function(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	tissue_names_arr <- c()
	per_gene_h2_arr <- c()
	per_gene_h2_se_arr <- c()	


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		non_med_rss_file <- paste0(tgfm_h2_results_dir, "tgfm_rss_likelihood_parallel_style_heritability_", trait_name, "_cis_heritable_gene_standardize_expr_True__robust_pleiotropic_prior_precision_temp.txt")
		non_med_rss_df <- read.table(non_med_rss_file, header=TRUE)

		h2_mean = non_med_rss_df$gamma_beta_b/non_med_rss_df$gamma_beta_a
		#h2_se1 = sqrt(tissue_df$gamma_alpha_b/((tissue_df$gamma_alpha_a)**2))
		h2_se = 1.0/sqrt(non_med_rss_df$gamma_beta_a/((non_med_rss_df$gamma_beta_b)**2))
	

		trait_names_arr <- c(trait_names_arr, trait_name)
		tissue_names_arr <- c(tissue_names_arr, "Non-mediated")

		per_gene_h2_arr <- c(per_gene_h2_arr, h2_mean[1])
		per_gene_h2_se_arr <- c(per_gene_h2_se_arr, h2_se[1])	


		# Load in jacknife data for this trait
		rss_file <- paste0(tgfm_h2_results_dir, "tgfm_rss_likelihood_parallel_style_heritability_", trait_name, "_cis_heritable_gene_standardize_expr_True_robust_tissue_specific_prior_precision_temp.txt")
		rss_df <- read.table(rss_file, header=TRUE)

		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- tissue_names[tissue_iter]
			tissue_name_raw <- tissue_names_raw[tissue_iter]
			tissue_df <- rss_df[as.character(rss_df$tissue)==tissue_name_raw,]


			h2_mean = tissue_df$gamma_alpha_b/tissue_df$gamma_alpha_a
			#h2_se1 = sqrt(tissue_df$gamma_alpha_b/((tissue_df$gamma_alpha_a)**2))
			h2_se = 1.0/sqrt(tissue_df$gamma_alpha_a/((tissue_df$gamma_alpha_b)**2))
	

			trait_names_arr <- c(trait_names_arr, trait_name)
			tissue_names_arr <- c(tissue_names_arr, tissue_name)
			per_gene_h2_arr <- c(per_gene_h2_arr, h2_mean)
			per_gene_h2_se_arr <- c(per_gene_h2_se_arr, h2_se)	
		}
	}

	# Put all into compact df
	df <- data.frame(trait=trait_names_arr, tissue=tissue_names_arr, per_gene_h2=per_gene_h2_arr, per_gene_h2_se=per_gene_h2_se_arr)

	return(df)	
}


extract_rss_per_gene_h2_df <- function(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw) {
	# Initialize arrays to save data
	trait_names_arr <- c()
	tissue_names_arr <- c()
	per_gene_h2_arr <- c()
	per_gene_h2_se_arr <- c()	


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]

		# Load in jacknife data for this trait
		rss_file <- paste0(tgfm_h2_results_dir, "tgfm_rss_likelihood_parallel_style_heritability_", trait_name, "_cis_heritable_gene_standardize_expr_True_robust_tissue_specific_prior_precision_temp.txt")
		rss_df <- read.table(rss_file, header=TRUE)

		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- tissue_names[tissue_iter]
			tissue_name_raw <- tissue_names_raw[tissue_iter]
			tissue_df <- rss_df[as.character(rss_df$tissue)==tissue_name_raw,]


			h2_mean = tissue_df$gamma_alpha_b/tissue_df$gamma_alpha_a
			#h2_se1 = sqrt(tissue_df$gamma_alpha_b/((tissue_df$gamma_alpha_a)**2))
			h2_se = 1.0/sqrt(tissue_df$gamma_alpha_a/((tissue_df$gamma_alpha_b)**2))
	

			trait_names_arr <- c(trait_names_arr, trait_name)
			tissue_names_arr <- c(tissue_names_arr, tissue_name)
			per_gene_h2_arr <- c(per_gene_h2_arr, h2_mean)
			per_gene_h2_se_arr <- c(per_gene_h2_se_arr, h2_se)	
		}
	}

	# Put all into compact df
	df <- data.frame(trait=trait_names_arr, tissue=tissue_names_arr, per_gene_h2=per_gene_h2_arr, per_gene_h2_se=per_gene_h2_se_arr)


	return(df)	
}


make_rss_ldsc_per_gene_h2_mean_se_for_single_trait <- function(rss_df, ldsc_df, trait_name) {
	tissue_arr <- as.character(c(as.character(rss_df$tissue), as.character(ldsc_df$tissue)))
	h2_arr <- c(rss_df$per_gene_h2, ldsc_df$per_gene_h2)
	model_type_arr <- c(rep("RSS", length(rss_df$per_gene_h2)), rep("LDSC", length(ldsc_df$per_gene_h2)))
	h2_lb_arr <- c(rss_df$per_gene_h2, ldsc_df$per_gene_h2 - 1.96*ldsc_df$per_gene_h2_se)
	h2_ub_arr <- c(rss_df$per_gene_h2, ldsc_df$per_gene_h2 + 1.96*ldsc_df$per_gene_h2_se)


	df <- data.frame(tissue=tissue_arr, h2=h2_arr, h2_lb=h2_lb_arr, h2_ub=h2_ub_arr, model_type=model_type_arr)

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")


	p<- ggplot(df, aes(x=tissue, y=h2, fill=model_type)) + 
 		 geom_bar(stat="identity", color="black", position=position_dodge()) +
 		 figure_theme() +
  		 geom_errorbar(aes(ymin=h2_lb, ymax=h2_ub), width=.2, position=position_dodge(.9)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		 theme(legend.position="none")+
  		 labs(x="", y="Per gene h2", title=trait_name, fill="") +
  		 theme(legend.position="bottom")
  	return(p)
}


make_rss_ldsc_gene_h2_mean_se_for_single_trait <- function(rss_df, ldsc_df, trait_name, num_genes_per_tissue) {
	tissue_order = as.character(rss_df$tissue)
	tissue_arr <- as.character(c(as.character(rss_df$tissue), as.character(ldsc_df$tissue)))
	h2_arr <- c(rss_df$per_gene_h2*num_genes_per_tissue, ldsc_df$per_gene_h2*num_genes_per_tissue)
	model_type_arr <- c(rep("RSS", length(rss_df$per_gene_h2)), rep("LDSC", length(ldsc_df$per_gene_h2)))
	h2_lb_arr <- c(rss_df$per_gene_h2, num_genes_per_tissue*ldsc_df$per_gene_h2 - 1.96*ldsc_df$per_gene_h2_se*num_genes_per_tissue)
	h2_ub_arr <- c(rss_df$per_gene_h2, num_genes_per_tissue*ldsc_df$per_gene_h2 + 1.96*ldsc_df$per_gene_h2_se*num_genes_per_tissue)


	total_rss_h2 = sum(rss_df$per_gene_h2*num_genes_per_tissue)
	total_ldsc_h2 = sum(ldsc_df$per_gene_h2*num_genes_per_tissue)

	df <- data.frame(tissue=tissue_arr, h2=h2_arr, h2_lb=h2_lb_arr, h2_ub=h2_ub_arr, model_type=model_type_arr)
	df$tissue = factor(df$tissue,levels=tissue_order)

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")


	p<- ggplot(df, aes(x=tissue, y=h2, fill=model_type)) + 
 		 geom_bar(stat="identity", color="black", position=position_dodge()) +
 		 figure_theme() +
  		 geom_errorbar(aes(ymin=h2_lb, ymax=h2_ub), width=.2, position=position_dodge(.9)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		 theme(legend.position="none")+
  		 labs(x="", y="Gene h2", title=paste0(trait_name, " / total rss h2: ", total_rss_h2, " / total ldsc h2: ", total_ldsc_h2), fill="") +
  		 theme(legend.position="bottom")
  	return(p)
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

make_fraction_med_per_trait_barplot_with_standard_errors <- function(df) {
	ordered_trait_names <- as.character(df$trait)
	ord <- order(df$fraction_h2_med)

	df$trait = factor(df$trait, levels=ordered_trait_names[ord])

	df$lb = df$fraction_h2_med - 1.96*df$fraction_h2_med_se
	df$ub = df$fraction_h2_med + 1.96*df$fraction_h2_med_se


	p<-ggplot(data=df, aes(x=fraction_h2_med, y=trait)) +
  		geom_bar(stat="identity") +
  		geom_errorbar(aes(xmin=lb, xmax=ub), width=.2, position=position_dodge(.9))  +
  		figure_theme()

  	return(p)
}


make_fraction_med_per_trait_barplot <- function(df) {
	ordered_trait_names <- as.character(df$trait)
	ord <- order(df$fraction_h2_med)

	df$trait = factor(df$trait, levels=ordered_trait_names[ord])



	p<-ggplot(data=df, aes(x=fraction_h2_med, y=trait)) +
  		geom_bar(stat="identity") +
  		figure_theme()

  	return(p)
}

scatterplot_of_h2_estimates <- function(x_val, y_val, x_axis_title, y_axis_title) {
	df <- data.frame(x_val=x_val, y_val=y_val)

	# Basic scatter plot
	p <- ggplot(df, aes(x=x_val, y=y_val)) + geom_point() +
	figure_theme() +
	labs(x=x_axis_title, y=y_axis_title) +
	geom_abline() 
	return(p)

}


######################
# Command line args
#######################
tissue_names_file = args[1]
trait_names_file = args[2]
num_genes_and_variants_dir = args[3]
tgfm_h2_results_dir = args[4]
standard_sldsc_results_dir = args[5]
standard_sldsc_processed_data_dir = args[6]
ldsc_baseline_ld_hg19_annotation_dir = args[7]
viz_dir = args[8]

# Load in gtex tissues
tissue_df <- read.table(tissue_names_file,header=TRUE)
tissue_names_raw <- as.character(tissue_df$pseudotissue_name)
tissue_names <- str_replace_all(tissue_names_raw, "-", "_")
num_genes_per_tissue <- get_number_of_genes_per_tissue(num_genes_and_variants_dir)
num_variants <- get_numer_of_variants(num_genes_and_variants_dir)


trait_names_df <- read.table(trait_names_file, header=TRUE,sep="\t")
# Get list of trait names
trait_names_df <- read.table(trait_names_file, header=TRUE,sep="\t")
trait_names <- as.character(trait_names_df$study_name)

if (FALSE) {
# Extract fraction of trait mediated heritability data based on standard sldsc
standard_sldsc_5_50_fraction_h2_med_expr_df = extract_standard_sldsc_fraction_of_h2_med_expression_df(standard_sldsc_results_dir, standard_sldsc_processed_data_dir, trait_names, tissue_names, tissue_names_raw, "l2.M_5_50")
standard_sldsc_all_fraction_h2_med_expr_df = extract_standard_sldsc_fraction_of_h2_med_expression_df(standard_sldsc_results_dir, standard_sldsc_processed_data_dir, trait_names, tissue_names, tissue_names_raw, "l2.M")
standard_sldsc_all_anno_5_50_fraction_h2_med_expr_df = extract_standard_sldsc_all_anno_fraction_of_h2_med_expression_df(standard_sldsc_results_dir, ldsc_baseline_ld_hg19_annotation_dir, trait_names, tissue_names, tissue_names_raw, "l2.M_5_50")
standard_sldsc_non_eqtl_cpp_removed_5_50_fraction_h2_med_expr_df = extract_standard_sldsc_non_eqtl_cpp_removed_fraction_of_h2_med_expression_df(standard_sldsc_results_dir, standard_sldsc_processed_data_dir, trait_names, tissue_names, tissue_names_raw, "l2.M_5_50")


# Extract fraction of trait mediated heritability data
rss_fraction_h2_med_expr_df <- extract_rss_fraction_of_h2_med_expression_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw, c(num_variants,num_genes_per_tissue))
ldsc_fraction_h2_med_expr_df <- extract_ldsc_fraction_of_h2_med_expression_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw, c(num_variants,num_genes_per_tissue))
ldsc_fraction_h2_med_expr_top_genes_df <- extract_ldsc_fraction_of_h2_med_expression_df2(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw, c(num_variants,num_genes_per_tissue))


# Fraction mediated bar plot
output_file <- paste0(viz_dir, "standard_sldsc_anno_non_eqtl_cpp_removed_5_50_fraction_mediated_per_trait.pdf")
fraction_med_per_trait_barplot = make_fraction_med_per_trait_barplot(standard_sldsc_non_eqtl_cpp_removed_5_50_fraction_h2_med_expr_df)
ggsave(fraction_med_per_trait_barplot, file=output_file, width=7.2, height=7.0, units="in")


# Fraction mediated
output_file <- paste0(viz_dir, "standard_sldsc_5_50_vs_tgfm_sldsc_fraction_mediated_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(standard_sldsc_5_50_fraction_h2_med_expr_df$fraction_h2_med, ldsc_fraction_h2_med_expr_df$fraction_h2_med, "Standard S-LDSC (5-50)\nh2_expr/(h2_expr+h2_geno)", "TGFM-LDSC\nh2_expr/(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "standard_sldsc_all_vs_tgfm_sldsc_fraction_mediated_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(standard_sldsc_all_fraction_h2_med_expr_df$fraction_h2_med, ldsc_fraction_h2_med_expr_df$fraction_h2_med, "Standard S-LDSC\nh2_expr/(h2_expr+h2_geno)", "TGFM-LDSC\nh2_expr/(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "tgfm_rss_vs_tgfm_sldsc_fraction_mediated_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(rss_fraction_h2_med_expr_df$fraction_h2_med, ldsc_fraction_h2_med_expr_df$fraction_h2_med, "TGFM-RSS\nh2_expr/(h2_expr+h2_geno)", "TGFM-LDSC\nh2_expr/(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "tgfm_sldsc_top_window_vs_tgfm_sldsc_fraction_mediated_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(ldsc_fraction_h2_med_expr_top_genes_df$fraction_h2_med, ldsc_fraction_h2_med_expr_df$fraction_h2_med, "TGFM-LDSC (Top window)\nh2_expr/(h2_expr+h2_geno)", "TGFM-LDSC\nh2_expr/(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "tgfm_sldsc_top_window_vs_tgfm_rss_fraction_mediated_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(ldsc_fraction_h2_med_expr_top_genes_df$fraction_h2_med, rss_fraction_h2_med_expr_df$fraction_h2_med, "TGFM-LDSC (Top window)\nh2_expr/(h2_expr+h2_geno)", "TGFM-RSS\nh2_expr/(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

# Total h2
output_file <- paste0(viz_dir, "standard_sldsc_5_50_vs_tgfm_sldsc_total_h2_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(standard_sldsc_5_50_fraction_h2_med_expr_df$total_h2, ldsc_fraction_h2_med_expr_df$total_h2, "Standard S-LDSC (5-50)\n(h2_expr+h2_geno)", "TGFM-LDSC\n(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "standard_sldsc_all_vs_tgfm_sldsc_total_h2_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(standard_sldsc_all_fraction_h2_med_expr_df$total_h2, ldsc_fraction_h2_med_expr_df$total_h2, "Standard S-LDSC\n(h2_expr+h2_geno)", "TGFM-LDSC\n(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "tgfm_rss_vs_tgfm_sldsc_total_h2_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(rss_fraction_h2_med_expr_df$total_h2, ldsc_fraction_h2_med_expr_df$total_h2, "TGFM-RSS\n(h2_expr+h2_geno)", "TGFM-LDSC\n(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "tgfm_sldsc_top_window_vs_tgfm_sldsc_total_h2_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(ldsc_fraction_h2_med_expr_top_genes_df$total_h2, ldsc_fraction_h2_med_expr_df$total_h2, "TGFM-SLDSC (Top window)\n(h2_expr+h2_geno)", "TGFM-LDSC\n(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "tgfm_sldsc_top_window_vs_tgfm_rss_total_h2_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(ldsc_fraction_h2_med_expr_top_genes_df$total_h2, rss_fraction_h2_med_expr_df$total_h2, "TGFM-SLDSC (Top window)\n(h2_expr+h2_geno)", "TGFM-RSS\n(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

# Fraction mediated vs total h2
output_file <- paste0(viz_dir, "tgfm_sldsc_total_h2_vs_tgfm_sldsc_fraction_mediated_scatter.pdf")
scatter <- scatterplot_of_h2_estimates(ldsc_fraction_h2_med_expr_df$total_h2, ldsc_fraction_h2_med_expr_df$fraction_h2_med, "TGFM-LDSC\n(h2_expr+h2_geno)", "TGFM-LDSC\nh2_expr/(h2_expr+h2_geno)")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")


# Fraction mediated bar plot
output_file <- paste0(viz_dir, "tgfm_sldsc_fraction_mediated_per_trait.pdf")
fraction_med_per_trait_barplot = make_fraction_med_per_trait_barplot_with_standard_errors(ldsc_fraction_h2_med_expr_df)
ggsave(fraction_med_per_trait_barplot, file=output_file, width=7.2, height=7.0, units="in")



# Fraction mediated bar plot
output_file <- paste0(viz_dir, "standard_sldsc_5_50_fraction_mediated_per_trait.pdf")
fraction_med_per_trait_barplot = make_fraction_med_per_trait_barplot(standard_sldsc_5_50_fraction_h2_med_expr_df)
ggsave(fraction_med_per_trait_barplot, file=output_file, width=7.2, height=7.0, units="in")

# Fraction mediated bar plot
output_file <- paste0(viz_dir, "standard_sldsc_all_anno_5_50_fraction_mediated_per_trait.pdf")
fraction_med_per_trait_barplot = make_fraction_med_per_trait_barplot(standard_sldsc_all_anno_5_50_fraction_h2_med_expr_df)
ggsave(fraction_med_per_trait_barplot, file=output_file, width=7.2, height=7.0, units="in")

}




#ordered_trait_names_gt = as.character(rss_fraction_h2_med_expr_df$trait)[rss_fraction_h2_med_expr_df$fraction_h2_med > .05]
ordered_trait_names_gt = trait_names

rss_per_gene_h2_df <- extract_rss_per_gene_h2_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw)
output_file <- paste0(viz_dir, "rss_fraction_h2_med_per_tissue_heatmap.pdf")
heatmap <- make_rss_variance_proportion_heatmap(rss_per_gene_h2_df, ordered_trait_names_gt)
ggsave(heatmap, file=output_file, width=7.2, height=7.0, units="in")

ldsc_per_gene_h2_df <- extract_ldsc_per_gene_h2_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw)
output_file <- paste0(viz_dir, "ldsc_fraction_h2_med_per_tissue_heatmap.pdf")
heatmap <- make_rss_variance_proportion_heatmap(ldsc_per_gene_h2_df, ordered_trait_names_gt)
ggsave(heatmap, file=output_file, width=7.2, height=7.0, units="in")



# Make heatmap showing z-scores of per gene trait heritabilities across tissues
output_file <- paste0(viz_dir , "tissue_z_score_heatmap_all_traits.pdf")
heatmap <- make_z_score_heatmap(tgfm_h2_results_dir, trait_names, tissue_names_raw)
ggsave(heatmap, file=output_file, width=13.2, height=7.0, units="in")

if (FALSE) {


output_file <- paste0(viz_dir, "fraction_mediated_per_trait.pdf")
fraction_med_per_trait_barplot = make_fraction_med_per_trait_barplot(rss_fraction_h2_med_expr_df)
ggsave(fraction_med_per_trait_barplot, file=output_file, width=7.2, height=7.0, units="in")


# Extract data
rss_per_gene_h2_plus_non_med_df <- extract_rss_per_gene_h2_plus_non_med_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw)
ldsc_per_gene_h2_plus_non_med_df <- extract_ldsc_per_gene_h2_plus_non_mediated_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw)

# Make gene h2 mean-se plot for single trait for both rss and ldsc
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(viz_dir, "rss_ldsc_gene_and_variants_h2_mean_se_", trait_name, ".pdf")
	rss_ldsc_per_gene_h2_mean_se = make_rss_ldsc_gene_h2_mean_se_for_single_trait(rss_per_gene_h2_plus_non_med_df[rss_per_gene_h2_plus_non_med_df$trait==trait_name,], ldsc_per_gene_h2_plus_non_med_df[ldsc_per_gene_h2_plus_non_med_df$trait==trait_name,], trait_name, c(num_variants,num_genes_per_tissue))
	ggsave(rss_ldsc_per_gene_h2_mean_se, file=output_file, width=7.2, height=4.0, units="in")
}



# Extract data
rss_per_gene_h2_df <- extract_rss_per_gene_h2_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw)
ldsc_per_gene_h2_df <- extract_ldsc_per_gene_h2_df(tgfm_h2_results_dir, trait_names, tissue_names, tissue_names_raw)


# Make gene h2 mean-se plot for single trait for both rss and ldsc
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(viz_dir, "rss_ldsc_gene_h2_mean_se_", trait_name, ".pdf")
	rss_ldsc_per_gene_h2_mean_se = make_rss_ldsc_gene_h2_mean_se_for_single_trait(rss_per_gene_h2_df[rss_per_gene_h2_df$trait==trait_name,], ldsc_per_gene_h2_df[ldsc_per_gene_h2_df$trait==trait_name,], trait_name, num_genes_per_tissue)
	ggsave(rss_ldsc_per_gene_h2_mean_se, file=output_file, width=7.2, height=4.0, units="in")
}


# Make per gene h2 mean-se plot for single trait for both rss and ldsc
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(viz_dir, "rss_ldsc_per_gene_h2_mean_se_", trait_name, ".pdf")
	rss_ldsc_per_gene_h2_mean_se = make_rss_ldsc_per_gene_h2_mean_se_for_single_trait(rss_per_gene_h2_df[rss_per_gene_h2_df$trait==trait_name,], ldsc_per_gene_h2_df[ldsc_per_gene_h2_df$trait==trait_name,], trait_name)
	ggsave(rss_ldsc_per_gene_h2_mean_se, file=output_file, width=7.2, height=4.0, units="in")
}

}

if (FALSE) {
#############
# Analysis specifically looking at LDSC style results
#############


# trait names
trait_df <- read.table(trait_names_file, header=TRUE)
#trait_df = trait_df[5:(dim(trait_df)[1]),]
trait_df = trait_df[as.character(trait_df$study_name) != "blood_WHITE_COUNT", ]
trait_df = trait_df[as.character(trait_df$study_name) != "biochemistry_HDLcholesterol", ]


trait_names_big <- as.character(trait_df$study_name)



# Make heatmap showing relative proportions of per gene trait heritabilities across tissues
output_file <- paste0(viz_dir , "tissue_variance_proportion_heatmap_all_traits.pdf")
heatmap <- make_variance_proportion_heatmap(tgfm_h2_results_dir, trait_names_big, tissue_names)
ggsave(heatmap, file=output_file, width=13.2, height=7.0, units="in")


# Make heatmap showing z-scores of per gene trait heritabilities across tissues
output_file <- paste0(viz_dir , "tissue_z_score_heatmap_all_traits.pdf")
heatmap <- make_z_score_heatmap(tgfm_h2_results_dir, trait_names_big, tissue_names_raw)
ggsave(heatmap, file=output_file, width=13.2, height=7.0, units="in")


trait_names_small <- rev(c("biochemistry_Cholesterol", "blood_EOSINOPHIL_COUNT", "blood_RBC_DISTRIB_WIDTH", "blood_RED_COUNT", "blood_LYMPHOCYTE_COUNT", "body_WHRadjBMIz", "bp_DIASTOLICadjMEDz", "cov_EDU_COLLEGE", "disease_ALLERGY_ECZEMA_DIAGNOSED"))
# Make heatmap showing relative proportions of per gene trait heritabilities across tissues
output_file <- paste0(viz_dir , "tissue_variance_proportion_heatmap.pdf")
heatmap <- make_variance_proportion_heatmap(tgfm_h2_results_dir, trait_names_small, tissue_names)
ggsave(heatmap, file=output_file, width=13.2, height=7.0, units="in")


# Make heatmap showing z-scores of per gene trait heritabilities across tissues
output_file <- paste0(viz_dir , "tissue_z_score_heatmap.pdf")
heatmap <- make_z_score_heatmap(tgfm_h2_results_dir, trait_names_small, tissue_names_raw)
ggsave(heatmap, file=output_file, width=13.2, height=7.0, units="in")



for (trait_iter in 1:length(trait_names_big)) {
	
	trait_name <- trait_names_big[trait_iter]
	se_plot <- jacknife_mean_standard_error_plot_for_single_trait(trait_name, tgfm_h2_results_dir, tissue_names_raw)
	output_file <- paste0(viz_dir , "gene_h2_with_standard_error_barplot_", trait_name, ".pdf")

	ggsave(se_plot, file=output_file, width=7.2, height=3.7, units="in")

}
}

