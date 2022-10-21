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
		variance_file = paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_learn_intercept_mean_estimates.txt")
		# load in df
		variance_df <- read.table(variance_file, header=TRUE)
		# filter to genes
		variance_df = variance_df[3:(dim(variance_df)[1]),]

		trait_tissue_names <- as.character(variance_df$Class_name)
		trait_var = variance_df$h2
		
		jacknife_file <- paste0(tgfm_h2_results_dir, "tgfm_ldsc_style_heritability_", trait_name, "_learn_intercept_jacknifed_mean_estimates.txt")
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





######################
# Command line args
#######################
tissue_names_file = args[1]
trait_names_file = args[2]
tgfm_h2_results_dir = args[3]
viz_dir = args[4]

# Load in gtex tissues
tissue_df <- read.table(tissue_names_file,header=TRUE)
tissue_names_raw <- as.character(tissue_df$pseudotissue_name)
tissue_names <- str_replace_all(tissue_names_raw, "-", "_")

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


