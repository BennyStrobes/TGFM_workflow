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

extract_fraction_of_h2_mediated_by_gene_expression_for_sparse_models <- function(full_trait_names, hr_trait_names, tgfm_sldsc_results_dir, methods, anno_size) {
	trait_arr <- c()
	method_arr <- c()
	med_h2_arr <- c()
	med_h2_se_arr <- c()
	med_h2_lb_arr <- c()
	med_h2_ub_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]
		hr_trait_name <- hr_trait_names[trait_iter]

		for (method_iter in 1:length(methods)) {
			method <- methods[method_iter]
			per_ele_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baselineLD_no_qtl_gene_adj_ld_scores_organized_", method, "_res.txt")
			tmp_df <- read.table(per_ele_h2_file, header=TRUE)

			partition_h2 = tmp_df$tau*anno_size
			geno_h2 = sum(partition_h2[1:93])
			expr_h2 = sum(partition_h2[94:length(partition_h2)])
	
			med_h2 = expr_h2/(expr_h2 + geno_h2)

			trait_arr <- c(trait_arr, hr_trait_name)
			method_arr <- c(method_arr, method)
			med_h2_arr <- c(med_h2_arr, med_h2)
			med_h2_se_arr <- c(med_h2_se_arr, 1.0)
			med_h2_lb_arr <- c(med_h2_lb_arr, med_h2)
			med_h2_ub_arr <- c(med_h2_ub_arr, med_h2)
		}
	}
	df <- data.frame(trait=trait_arr, method=method_arr, fraction_h2=med_h2_arr, med_h2_se=med_h2_se_arr, fraction_h2_lb=med_h2_lb_arr, fraction_h2_ub=med_h2_ub_arr)
	return(df)	
}


extract_fraction_of_h2_mediated_by_gene_expression_for_several_methods <- function(full_trait_names, hr_trait_names, tgfm_sldsc_results_dir, methods) {
	trait_arr <- c()
	method_arr <- c()
	med_h2_arr <- c()
	med_h2_se_arr <- c()
	med_h2_lb_arr <- c()
	med_h2_ub_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]
		hr_trait_name <- hr_trait_names[trait_iter]

		for (method_iter in 1:length(methods)) {
			method <- methods[method_iter]
			expr_med_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", method, "_gene_adj_ld_scores_h2_5_50_med.txt")
			tmp_df <- read.table(expr_med_h2_file, header=TRUE)

			trait_arr <- c(trait_arr, hr_trait_name)
			method_arr <- c(method_arr, method)
			med_h2_arr <- c(med_h2_arr, tmp_df$h2_med[1])
			med_h2_se_arr <- c(med_h2_se_arr, tmp_df$h2_med_se[1])
			med_h2_lb_arr <- c(med_h2_lb_arr, tmp_df$h2_med[1] - 1.96*tmp_df$h2_med_se[1])
			med_h2_ub_arr <- c(med_h2_ub_arr, tmp_df$h2_med[1] + 1.96*tmp_df$h2_med_se[1])
		}
	}
	df <- data.frame(trait=trait_arr, method=method_arr, fraction_h2=med_h2_arr, med_h2_se=med_h2_se_arr, fraction_h2_lb=med_h2_lb_arr, fraction_h2_ub=med_h2_ub_arr)
	return(df)
}

make_fraction_mediated_h2_barplot <- function(df, ordered_levels) {
	sub_df <- df[as.character(df$method)=="baselineLD_no_qtl",]
	initial_trait_ordering = as.character(sub_df$trait)
	ord <- order(sub_df$fraction_h2)
	new_trait_ordering = initial_trait_ordering[ord]

	df$trait = factor(df$trait, levels=new_trait_ordering)
	df$method = factor(df$method, levels=ordered_levels)

	p<-ggplot(data=df, aes(x=trait, y=fraction_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_h2_lb, ymax=fraction_h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11)) +
  		labs(title="", x="", y="fraction h2 mediated\nby genetic gene expression", fill="") +
  		theme(legend.position="top")

  	return(p)	

}

load_in_ldsc_results <- function(tgfm_sldsc_results_dir, full_trait_names, hr_trait_names, methods) {
	trait_arr <- c()
	method_arr <- c()
	tissue_arr <- c()
	h2_arr <- c()
	h2_se_arr <- c()
	h2_z_arr <- c()
	h2_lb_arr <- c()
	h2_ub_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]
		hr_trait_name <- hr_trait_names[trait_iter]

		for (method_iter in 1:length(methods)) {
			method <- methods[method_iter]
			per_ele_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", method, "_gene_adj_ld_scores_organized_res.txt")
			per_ele_h2_df <- read.table(per_ele_h2_file, header=TRUE)
			num_anno = length(per_ele_h2_df$Annotation)
			start_index = (1:num_anno)[per_ele_h2_df$Annotation == "Adipose_Subcutaneous"]
			
			per_ele_h2_df = per_ele_h2_df[start_index:num_anno,]

			n_tiss = length(per_ele_h2_df$Annotation)
			tissue_arr <- c(tissue_arr, as.character(per_ele_h2_df$Annotation))
			trait_arr <- c(trait_arr, rep(hr_trait_name, n_tiss))
			method_arr <- c(method_arr, rep(method, n_tiss))
			h2_arr <- c(h2_arr, per_ele_h2_df$tau)
			h2_se_arr <- c(h2_se_arr, per_ele_h2_df$tau_se)
			h2_z_arr <- c(h2_z_arr, per_ele_h2_df$tau_z)
			h2_lb_arr <- c(h2_lb_arr, per_ele_h2_df$tau - 1.96*per_ele_h2_df$tau_se)
			h2_ub_arr <- c(h2_ub_arr, per_ele_h2_df$tau + 1.96*per_ele_h2_df$tau_se)
		}
	}

	df <- data.frame(trait=trait_arr, method=method_arr, tissue=tissue_arr, h2=h2_arr, h2_se=h2_se_arr, h2_z=h2_z_arr, h2_lb=h2_lb_arr, h2_ub=h2_ub_arr)
	return(df)
}

load_in_rss_results <- function(tgfm_heritability_results_dir, full_trait_names, hr_trait_names, methods) {
	trait_arr <- c()
	method_arr <- c()
	tissue_arr <- c()
	h2_arr <- c()
	h2_se_arr <- c()
	h2_z_arr <- c()
	h2_lb_arr <- c()
	h2_ub_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]
		hr_trait_name <- hr_trait_names[trait_iter]

		for (method_iter in 1:length(methods)) {
			method <- methods[method_iter]
			per_ele_h2_file <- paste0(tgfm_heritability_results_dir, "tgfm_rss_likelihood_parallel_style_heritability_", trait_name, "_cis_heritable_gene_standardize_expr_True_robust_tissue_specific_prior_precision_temp.txt")
			per_ele_h2_df <- read.table(per_ele_h2_file, header=TRUE)
			
			tissue_arr <- c(tissue_arr, as.character(per_ele_h2_df$tissue))
			n_tiss = length(per_ele_h2_df$tissue)
			trait_arr <- c(trait_arr, rep(hr_trait_name, n_tiss))
			method_arr <- c(method_arr, rep(method, n_tiss))
			h2_arr <- c(h2_arr, 1.0/per_ele_h2_df$expected_precision)
			h2_se_arr <- c(h2_se_arr, 0.0)
			h2_z_arr <- c(h2_z_arr, 0.0)
			h2_lb_arr <- c(h2_lb_arr, 1.0/per_ele_h2_df$expected_precision)
			h2_ub_arr <- c(h2_ub_arr, 1.0/per_ele_h2_df$expected_precision)
		}
	}

	df <- data.frame(trait=trait_arr, method=method_arr, tissue=tissue_arr, h2=h2_arr, h2_se=h2_se_arr, h2_z=h2_z_arr, h2_lb=h2_lb_arr, h2_ub=h2_ub_arr)
	return(df)
}

load_in_sparse_ldsc_results <- function(tgfm_sldsc_results_dir, full_trait_names, hr_trait_names, methods) {
	trait_arr <- c()
	method_arr <- c()
	tissue_arr <- c()
	h2_arr <- c()
	h2_se_arr <- c()
	h2_z_arr <- c()
	h2_lb_arr <- c()
	h2_ub_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]
		hr_trait_name <- hr_trait_names[trait_iter]

		for (method_iter in 1:length(methods)) {
			method <- methods[method_iter]
			per_ele_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baselineLD_no_qtl_gene_adj_ld_scores_organized_", method, "_res.txt")
			per_ele_h2_df <- read.table(per_ele_h2_file, header=TRUE)
			num_anno = length(per_ele_h2_df$Annotation)
			start_index = (1:num_anno)[per_ele_h2_df$Annotation == "Adipose_Subcutaneous"]
			
			per_ele_h2_df = per_ele_h2_df[start_index:num_anno,]

			n_tiss = length(per_ele_h2_df$Annotation)
			tissue_arr <- c(tissue_arr, as.character(per_ele_h2_df$Annotation))
			trait_arr <- c(trait_arr, rep(hr_trait_name, n_tiss))
			method_arr <- c(method_arr, rep(method, n_tiss))
			h2_arr <- c(h2_arr, per_ele_h2_df$tau)
			h2_se_arr <- c(h2_se_arr, per_ele_h2_df$tau_se)
			h2_z_arr <- c(h2_z_arr, per_ele_h2_df$tau_z)
			h2_lb_arr <- c(h2_lb_arr, per_ele_h2_df$tau - 1.96*per_ele_h2_df$tau_se)
			h2_ub_arr <- c(h2_ub_arr, per_ele_h2_df$tau + 1.96*per_ele_h2_df$tau_se)
		}
	}

	df <- data.frame(trait=trait_arr, method=method_arr, tissue=tissue_arr, h2=h2_arr, h2_se=h2_se_arr, h2_z=h2_z_arr, h2_lb=h2_lb_arr, h2_ub=h2_ub_arr)
	return(df)
}

make_per_predictor_se_barplot <- function(df, trait_name, levels_arr) {


	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")

	df$tissue = str_replace_all(as.character(df$tissue), "Expression_", "")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	df$method = factor(df$method, levels=levels_arr)

	#print(df)

	p<-ggplot(data=df, aes(x=tissue, y=h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=h2_lb, ymax=h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=trait_name, x="") +
  		theme(legend.position="top")

  	return(p)
}

make_coef_score_heatmap <- function(df, ordered_trait_names) {
	df$tissue = str_replace_all(as.character(df$tissue), "Expression_", "")
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$trait = factor(df$trait, ordered_trait_names)
	#df$h2_z = abs(df$h2_z)
	df$h2[df$h2 < 0.0] = 0.0
	p <- ggplot(df, aes(x = tissue, y = trait, fill = h2)) +
  		geom_tile() +
  		theme(text = element_text(size=11), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="Expression-mediated\ntrait heritability per gene",x="Tissue", y="GWAS trait",title="")

  	return(p)
}

make_z_score_heatmap <- function(df, ordered_trait_names, threshold=2.0) {
	df$tissue = str_replace_all(as.character(df$tissue), "Expression_", "")
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$trait = factor(df$trait, ordered_trait_names)
	#df$h2_z = abs(df$h2_z)
	df$h2_z[df$h2_z < 0.0] = 0.0
	df$h2_z[df$h2_z < threshold] = 0.0
	p <- ggplot(df, aes(x = tissue, y = trait, fill = h2_z)) +
  		geom_tile() +
  		theme(text = element_text(size=11), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="Expression-mediated\ntrait heritability z-score",x="Tissue", y="GWAS trait",title="")

  	return(p)
}

 extract_anno_size <- function(processed_tgfm_sldsc_data_dir) {
 	for (chrom_num in 1:22) {
 		chrom_filer = paste0(processed_tgfm_sldsc_data_dir, "baselineLD_no_qtl_gene_adj_ld_scores.", chrom_num, ".l2.M_5_50")
 		tmp_data = read.table(chrom_filer,header=FALSE)
 		if (chrom_num == 1) {
 			counter = as.numeric(tmp_data[1,])
 		} else {
 			counter = counter + as.numeric(tmp_data[1,])
 		}
 	}
 	return(counter)
}

std_mean <- function(x) sd(x)/sqrt(length(x))



############
# command line args
tissue_names_file = args[1]
trait_file = args[2]
processed_tgfm_sldsc_data_dir = args[3]
tgfm_sldsc_results_dir = args[4]
viz_dir = args[5]
tgfm_heritability_results_dir = args[6]


print(trait_file)
# Get list of independent traits
full_trait_names <- c("body_WHRadjBMIz", "body_BMIz", "body_HEIGHTz", "bmd_HEEL_TSCOREz", "blood_MEAN_PLATELET_VOL", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "biochemistry_LDLdirect", "bp_DIASTOLICadjMEDz", "disease_CARDIOVASCULAR", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "cov_EDU_COLLEGE", "biochemistry_VitaminD")
hr_trait_names <- c("WHR-adj_BMI","BMI", "Height", "BMD", "Platelet_vol", "Mean_corp_hemoglobin", "Monocyte_count", "Reticulocyte_count", "Pigment_hair", "FEV1FVCz", "Balding", "Cholesterol", "LDLdirect", "Diastolic_bp", "Cardiovascular_disease", "FVC", "Menarche_age", "Eczema", "College_education", "VitaminD")


#full_trait_names <- c("body_WHRadjBMIz", "body_BMIz", "body_HEIGHTz", "bmd_HEEL_TSCOREz", "blood_MEAN_PLATELET_VOL", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "cov_EDU_COLLEGE", "biochemistry_VitaminD")
#hr_trait_names <- c("WHR-adj_BMI","BMI", "Height", "BMD", "Platelet_vol", "Mean_corp_hemoglobin", "Monocyte_count", "Reticulocyte_count", "Pigment_hair", "FEV1FVCz", "Balding", "Cholesterol", "Diastolic_bp", "FVC", "Menarche_age", "Eczema", "College_education", "VitaminD")

# Load in gtex tissues
tissue_df <- read.table(tissue_names_file,header=TRUE)
tissue_names_raw <- as.character(tissue_df$pseudotissue_name)
tissue_names <- str_replace_all(tissue_names_raw, "-", "_")


# Extract fraction of heritability mediated by gene expression in each tissue
methods <- c("genotype_intercept", "baseline_no_qtl", "baselineLD_no_qtl")
expr_med_frac_df = extract_fraction_of_h2_mediated_by_gene_expression_for_several_methods(full_trait_names, hr_trait_names, tgfm_sldsc_results_dir, methods)


# Make fraction mediated h2 barplot for two models
fraction_mediated_barplot <- make_fraction_mediated_h2_barplot(expr_med_frac_df, methods)
output_file <- paste0(viz_dir, "fraction_h2_mediated_se_barplot.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")

# Extract anno size vector (across annotations)
anno_size = extract_anno_size(processed_tgfm_sldsc_data_dir)
methods <- c("sparse_ard_no_geno_regularization")
expr_med_frac_sparse_model_frac_df = extract_fraction_of_h2_mediated_by_gene_expression_for_sparse_models(full_trait_names, hr_trait_names, tgfm_sldsc_results_dir, methods, anno_size)

# Make fraction mediated h2 barplot for sparse and non-sparse models
tmp <- expr_med_frac_df[expr_med_frac_df$method=="baselineLD_no_qtl",]
cat_df = bind_rows(tmp, expr_med_frac_sparse_model_frac_df)
fraction_mediated_barplot <- make_fraction_mediated_h2_barplot(cat_df, c("baselineLD_no_qtl", "sparse_ard_no_geno_regularization"))
output_file <- paste0(viz_dir, "fraction_h2_mediated_se_sparse_vs_nonsparse_barplot.pdf")
ggsave(fraction_mediated_barplot, file=output_file, width=7.2, height=5.5, units="in")




full_trait_names <- c("body_WHRadjBMIz", "body_BMIz", "body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "lung_FEV1FVCzSMOKE", "biochemistry_Cholesterol", "biochemistry_LDLdirect", "bp_DIASTOLICadjMEDz", "disease_CARDIOVASCULAR", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "cov_EDU_COLLEGE", "biochemistry_VitaminD")
hr_trait_names <- c("WHR-adj_BMI","BMI", "Height", "Platelet_vol", "Mean_corp_hemoglobin", "Monocyte_count", "Reticulocyte_count", "FEV1FVCz", "Cholesterol", "LDLdirect", "Diastolic_bp", "Cardiovascular_disease", "FVC", "Menarche_age", "Eczema", "College_education", "VitaminD")


methods <- c("genotype_intercept", "baseline_no_qtl", "baselineLD_no_qtl")
ldsc_df <- load_in_ldsc_results(tgfm_sldsc_results_dir, full_trait_names, hr_trait_names, methods)


# Make z-score boxplot based on both
subset_df = ldsc_df[ldsc_df$method=="baselineLD_no_qtl",]
z_score_heatmap <- make_z_score_heatmap(subset_df, hr_trait_names, threshold=2.0)
output_file <- paste0(viz_dir, "baselineLD_nonneg_tau_z_score_heatmap.pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=6.0, units="in")


subset_df = ldsc_df[ldsc_df$method=="genotype_intercept",]
z_score_heatmap <- make_z_score_heatmap(subset_df, hr_trait_names, threshold=2.0)
output_file <- paste0(viz_dir, "genotype_intercept_nonneg_tau_z_score_heatmap.pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=6.0, units="in")



for (trait_iter in 1:length(hr_trait_names)) {
	hr_trait_name <- hr_trait_names[trait_iter]
	# Make se bar plot show per predictor h2
	tmp_ldsc_df = ldsc_df[ldsc_df$trait==hr_trait_name, ]

	predictor_se_barplot <- make_per_predictor_se_barplot(tmp_ldsc_df, hr_trait_name, c("genotype_intercept", "baseline_no_qtl", "baselineLD_no_qtl"))
	output_file <- paste0(viz_dir, "trait_specific_", hr_trait_name, "_per_predictor_tau_se_barplot.pdf")
	ggsave(predictor_se_barplot, file=output_file, width=8.6, height=4.0, units="in")

}


methods <- c("rss")
rss_df <- load_in_rss_results(tgfm_heritability_results_dir, full_trait_names, hr_trait_names, methods)

methods <- c("sparse_ard_no_geno_regularization")
sparse_ldsc_df <- load_in_sparse_ldsc_results(tgfm_sldsc_results_dir, full_trait_names, hr_trait_names, methods)


z_score_heatmap <- make_coef_score_heatmap(sparse_ldsc_df, hr_trait_names)
output_file <- paste0(viz_dir, "sparse_baselineld_nonneg_tau_coef_score_heatmap.pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=6.0, units="in")

z_score_heatmap <- make_coef_score_heatmap(rss_df, hr_trait_names)
output_file <- paste0(viz_dir, "rss_tau_coef_score_heatmap.pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=6.0, units="in")


for (trait_iter in 1:length(hr_trait_names)) {
	hr_trait_name <- hr_trait_names[trait_iter]
	# Make se bar plot show per predictor h2
	tmp_sparse_ldsc_df = sparse_ldsc_df[sparse_ldsc_df$trait==hr_trait_name, ]
	tmp_ldsc_df = ldsc_df[ldsc_df$trait==hr_trait_name, ]
	tmp_ldsc_df = tmp_ldsc_df[tmp_ldsc_df$method=="baselineLD_no_qtl", ]

	cat_df = bind_rows(tmp_ldsc_df, tmp_sparse_ldsc_df)


	predictor_se_barplot <- make_per_predictor_se_barplot(cat_df, hr_trait_name, c("baselineLD_no_qtl", "sparse_ard_no_geno_regularization"))
	output_file <- paste0(viz_dir, "trait_specific_", hr_trait_name, "_per_predictor_sparse_vs_not_tau_se_barplot.pdf")
	ggsave(predictor_se_barplot, file=output_file, width=8.6, height=4.0, units="in")
}




methods <- c("sparse_ard_no_geno_regularization", "sparse_nonneg_ard_no_geno_regularization")
sparse_ldsc_df <- load_in_sparse_ldsc_results(tgfm_sldsc_results_dir, full_trait_names, hr_trait_names, methods)

for (trait_iter in 1:length(hr_trait_names)) {
	hr_trait_name <- hr_trait_names[trait_iter]
	# Make se bar plot show per predictor h2
	tmp_ldsc_df = sparse_ldsc_df[sparse_ldsc_df$trait==hr_trait_name, ]

	predictor_se_barplot <- make_per_predictor_se_barplot(tmp_ldsc_df, hr_trait_name, methods)
	output_file <- paste0(viz_dir, "trait_specific_", hr_trait_name, "_per_sparse_predictor_tau_se_barplot.pdf")
	ggsave(predictor_se_barplot, file=output_file, width=8.6, height=4.0, units="in")
}
