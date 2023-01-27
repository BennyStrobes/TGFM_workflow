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


load_in_tgfm_sldsc_coefficients <- function(tgfm_sldsc_results_dir, trait_name, models, regularization_weights) {
	coef_arr <- c()
	coef_se_arr <- c()
	tissue_arr <- c()
	trait_arr <- c()
	model_arr <- c()
	regularization_weight_arr <- c()
	for (model_iter in 1:length(models)) {
		for (reg_iter in 1:length(regularization_weights)) {
			model_name <- models[model_iter]
			reg_weight = regularization_weights[reg_iter]
			filer <- paste0(tgfm_sldsc_results_dir, trait_name, "_baselineLD_no_qtl_pmces_gene_adj_ld_scores_gene_weighted_organized_", reg_weight, "_sparse_", model_name, "_res.txt")
			
			df <- read.table(filer, header=TRUE)
			df_sub = df[94:dim(df)[1],]
			

			coef_arr <- c(coef_arr, df_sub$tau)
			coef_se_arr <- c(coef_se_arr, df_sub$tau_se)
			tissue_arr <- c(tissue_arr, as.character(df_sub$Annotation))
			nn = length(df_sub$tau)
			trait_arr <- c(trait_arr, rep(trait_name, nn))
			model_arr <- c(model_arr, rep(model_name, nn))
			regularization_weight_arr <- c(regularization_weight_arr, rep(reg_weight, nn))

		}
	}
	df <- data.frame(tau=coef_arr, tau_se=coef_se_arr, tissue=tissue_arr, trait=trait_arr, model=model_arr, regularization_weight=regularization_weight_arr)
	
	return(df)
}


per_tissue_barplot_for_range_of_regularizations <- function(df, model_name, regularization_weights, trait_name) {
	df$regularization_weight = factor(df$regularization_weight, levels=regularization_weights)
	df$tau[df$tau < 0.0] = 0.0
	p <- ggplot(data=df, aes(x=tissue, y=tau, fill=regularization_weight)) +
		geom_bar(stat="identity", color="black", position=position_dodge())+
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		labs(x="",y="TGFM-SLDSC tau", title=paste0(trait_name, " / ", model_name), fill="regularization weight") +
  		theme(legend.position="top")

  	return(p)
}

make_ld_score_annotation_correlation_heatmap <- function(preprocessed_tgfm_sldsc_data_dir) {
	# Data file names
	matrix_file <- paste0(preprocessed_tgfm_sldsc_data_dir, "baselineLD_no_qtl_pmces_gene_adj_ld_scores_ld_score_annotation_correlation_matrix_subset.txt")
	row_names_file <- paste0(preprocessed_tgfm_sldsc_data_dir, "baselineLD_no_qtl_pmces_gene_adj_ld_scores_ld_score_annotation_correlation_matrix_subset_row_names.txt")
	col_names_file <- paste0(preprocessed_tgfm_sldsc_data_dir, "baselineLD_no_qtl_pmces_gene_adj_ld_scores_ld_score_annotation_correlation_matrix_subset_col_names.txt")
	# Load in data
	row_names <- as.character(read.table(row_names_file, header=FALSE)$V1)
	col_names <- as.character(read.table(col_names_file, header=FALSE)$V1)
	corr_mat <- as.matrix(read.table(matrix_file, header=FALSE))
	rownames(corr_mat) = row_names
	colnames(corr_mat) = col_names


	row_names_ord <- hclust( dist(scale((corr_mat)), method = "euclidean"), method = "ward.D" )$order
	col_names_ord <- hclust( dist(scale(t(corr_mat)), method = "euclidean"), method = "ward.D" )$order


    melted_mat <- melt(corr_mat)
    #colnames(melted_mat) <- c("Covariate", "PC","PVE")
    colnames(melted_mat) <- c("BaselineLD_Annotation", "eQTL_tissue","R")
    melted_mat$BaselineLD_Annotation = factor(melted_mat$BaselineLD_Annotation, levels=row_names[row_names_ord])
    #melted_mat$eQTL_tissue = factor(melted_mat$eQTL_tissue, levels=col_names[col_names_ord])

    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=eQTL_tissue, y=BaselineLD_Annotation)) + geom_tile(aes(fill=R)) + scale_fill_gradient2(midpoint=-.0005, guide="colorbar")
    #heatmap <- heatmap + theme(text = element_text(size=8), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))
    heatmap <- heatmap + figure_theme() + theme(axis.text.x=element_text(angle=90, vjust=.5)) + theme(text = element_text(size=8)) + theme(axis.text=element_text(size=9))

   	return(heatmap)
}

load_in_per_tissue_tau_df_for_sparse_and_non_sparse_models <- function(trait_name, tissue_version, gene_version, annotation_version, sparse_model_name, regularization_weights, tgfm_sldsc_results_dir) {
	tissue_vec <- c()
	tau_vec <- c()
	tau_se_vec <- c()
	model_type_vec <- c()

	expr_med_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", annotation_version, "_", gene_version, "_", tissue_version, "_pmces_gene_adj_ld_scores_organized_res.txt")
	t_df <- read.table(expr_med_h2_file, header=TRUE, sep="\t")
	first_eqtl_row_num = which(as.character(t_df$Annotation)=="Adipose_Subcutaneous")
	eqtl_df <- t_df[first_eqtl_row_num:dim(t_df)[1],]
	n_tiss = dim(eqtl_df)[1]
	tissue_vec <- c(tissue_vec, as.character(eqtl_df$Annotation))
	tau_vec <- c(tau_vec, eqtl_df$tau)
	tau_se_vec <- c(tau_se_vec, eqtl_df$tau_se)
	model_type_vec <- c(model_type_vec, rep("tgfm-sldsc", n_tiss))


	for (regularization_weight_iter in 1:length(regularization_weights)) {
		regularization_weight = regularization_weights[regularization_weight_iter]
		sparse_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", annotation_version, "_", gene_version, "_", tissue_version, "_pmces_gene_adj_ld_scores_organized_", regularization_weight, "_sparse_", sparse_model_name, "_res.txt")
		t_df <- read.table(sparse_file, header=TRUE, sep="\t")
		first_eqtl_row_num = which(as.character(t_df$Annotation)=="Adipose_Subcutaneous")
		eqtl_df <- t_df[first_eqtl_row_num:dim(t_df)[1],]
		n_tiss = dim(eqtl_df)[1]
		tissue_vec <- c(tissue_vec, as.character(eqtl_df$Annotation))
		tau_vec <- c(tau_vec, eqtl_df$tau)
		tau_se_vec <- c(tau_se_vec, eqtl_df$tau_se)
		model_type_vec <- c(model_type_vec, rep(paste0(sparse_model_name,"_", regularization_weight), n_tiss))
	}



	df <- data.frame(tissue=as.character(tissue_vec), tau=tau_vec, tau_se=tau_se_vec, model_type=model_type_vec)

	return(df)
}

load_in_per_tissue_tau_df <- function(trait_names, tissue_versions, gene_versions, annotation_versions, tgfm_sldsc_results_dir) {
	trait_vec <- c()
	tissue_version_vec <- c()
	tissue_vec <- c()
	gene_vec <- c()
	annotation_vec <- c()
	tau_vec <- c()
	tau_se_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		for (tissue_version_iter in 1:length(tissue_versions)) {
			tissue_version = tissue_versions[tissue_version_iter]
			for (gene_version_iter in 1:length(gene_versions)) {
				gene_version <- gene_versions[gene_version_iter]
				for (annotation_version_iter in 1:length(annotation_versions)) {
					annotation_version = annotation_versions[annotation_version_iter]

					expr_med_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", annotation_version, "_", gene_version, "_", tissue_version, "_pmces_gene_adj_ld_scores_organized_res.txt")

					t_df <- read.table(expr_med_h2_file, header=TRUE, sep="\t")

					first_eqtl_row_num = which(as.character(t_df$Annotation)=="Adipose_Subcutaneous")
					eqtl_df <- t_df[first_eqtl_row_num:dim(t_df)[1],]
					n_tiss = dim(eqtl_df)[1]

					trait_vec <- c(trait_vec, rep(trait_name, n_tiss))
					tissue_version_vec <- c(tissue_version_vec, rep(tissue_version, n_tiss))
					tissue_vec <- c(tissue_vec, as.character(eqtl_df$Annotation))
					gene_vec <- c(gene_vec, rep(gene_version, n_tiss))
					annotation_vec <- c(annotation_vec, rep(annotation_version, n_tiss))
					tau_vec <- c(tau_vec, eqtl_df$tau)
					tau_se_vec <- c(tau_se_vec, eqtl_df$tau_se)
				}
			}
		}
	}

	df <- data.frame(trait=as.character(trait_vec), tissue=as.character(tissue_vec), tissue_version=as.character(tissue_version_vec), gene_version=as.character(gene_vec), annotation_model=as.character(annotation_vec), tau=tau_vec, tau_se=tau_se_vec)

	return(df)
}




load_in_total_expression_mediated_h2_df <- function(trait_names, tissue_versions, gene_versions, annotation_versions, tgfm_sldsc_results_dir) {
	trait_vec <- c()
	tissue_vec <- c()
	gene_vec <- c()
	annotation_vec <- c()
	expr_med_vec <- c()
	expr_med_se_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		for (tissue_version_iter in 1:length(tissue_versions)) {
			tissue_version = tissue_versions[tissue_version_iter]
			for (gene_version_iter in 1:length(gene_versions)) {
				gene_version <- gene_versions[gene_version_iter]
				for (annotation_version_iter in 1:length(annotation_versions)) {
					annotation_version = annotation_versions[annotation_version_iter]

					expr_med_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", annotation_version, "_", gene_version, "_", tissue_version, "_pmces_gene_adj_ld_scores_h2_5_50_med.txt")

					t_df <- read.table(expr_med_h2_file, header=TRUE, sep="\t")

					trait_vec <- c(trait_vec, trait_name)
					tissue_vec <- c(tissue_vec, tissue_version)
					gene_vec <- c(gene_vec, gene_version)
					annotation_vec <- c(annotation_vec, annotation_version)
					expr_med_vec <- c(expr_med_vec, t_df$h2_med[1])
					expr_med_se_vec <- c(expr_med_se_vec, t_df$h2_med_se[1])

				}
			}
		}
	}

	df <- data.frame(trait=as.character(trait_vec), tissue_version=as.character(tissue_vec), gene_version=as.character(gene_vec), annotation_model=as.character(annotation_vec), expression_mediated_h2=expr_med_vec, expression_mediated_h2_se=expr_med_se_vec)

	return(df)
}


make_fraction_mediated_h2_barplot_stratefied_by_gene_and_tissue_models <- function(df, titler) {
	frac_lb_arr <- c(df$expression_mediated_h2 - 1.96*df$expression_mediated_h2_se)
	frac_ub_arr <- c(df$expression_mediated_h2 + 1.96*df$expression_mediated_h2_se)

	new_version = paste0(df$gene_version,"_", df$tissue_version)

	df2 <- data.frame(trait=df$trait,fraction_h2=df$expression_mediated_h2, fraction_h2_lb=frac_lb_arr, fraction_h2_ub=frac_ub_arr, method=new_version)

	df2$trait <- recode(df2$trait, pigment_HAIR="Hair pigment", body_BALDING1="Balding1",bmd_HEEL_TSCOREz="Heel t-score", biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")
	#df2$method = factor(df2$method, levels=ordered_annotation_models)

	tmp_df = df2[as.character(df2$method)==unique(df2$method)[1],]
	ord <- order(tmp_df$fraction_h2)
	df2$trait <- factor(df2$trait, levels=as.character(tmp_df$trait)[ord])
	p<-ggplot(data=df2, aes(x=trait, y=fraction_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_h2_lb, ymax=fraction_h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=titler, x="", y="Fraction h2\nmediated by gene expression", fill="") +
  		theme(legend.position="bottom") +guides(fill=guide_legend(nrow=2,byrow=TRUE))

  	return(p)	
}

make_fraction_mediated_h2_barplot_stratefied_by_annotation_models <- function(df, titler, ordered_annotation_models) {
	frac_lb_arr <- c(df$expression_mediated_h2 - 1.96*df$expression_mediated_h2_se)
	frac_ub_arr <- c(df$expression_mediated_h2 + 1.96*df$expression_mediated_h2_se)


	df2 <- data.frame(trait=df$trait,fraction_h2=df$expression_mediated_h2, fraction_h2_lb=frac_lb_arr, fraction_h2_ub=frac_ub_arr, method=df$annotation_model)
	#df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")
	df2$trait <- recode(df2$trait, pigment_HAIR="Hair pigment", body_BALDING1="Balding1",bmd_HEEL_TSCOREz="Heel t-score", biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")
	df2$method = factor(df2$method, levels=ordered_annotation_models)

	tmp_df = df2[as.character(df2$method)==ordered_annotation_models[1],]
	ord <- order(tmp_df$fraction_h2)
	df2$trait <- factor(df2$trait, levels=as.character(tmp_df$trait)[ord])

	p<-ggplot(data=df2, aes(x=trait, y=fraction_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_h2_lb, ymax=fraction_h2_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=titler, x="", y="Fraction h2\nmediated by gene expression", fill="") +
  		theme(legend.position="bottom")

  	return(p)	

}

make_per_tissue_tau_barplot_stratefied_by_gene_models <- function(df, titler, ordered_annotation_models) {
	frac_lb_arr <- c(df$tau - 1.96*df$tau_se)
	frac_ub_arr <- c(df$tau + 1.96*df$tau_se)

	#new_version = paste0(df$gene_version,"_", df$tissue_version)

	df2 <- data.frame(tissue=df$tissue,tau=df$tau, tau_lb=frac_lb_arr, tau_ub=frac_ub_arr, method=factor(df$gene_version))


	#df2$method = factor(df2$method, levels=ordered_annotation_models)

	p<-ggplot(data=df2, aes(x=tissue, y=tau, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=tau_lb, ymax=tau_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=titler, x="", y="TGFM-SLDSC Tau", fill="") +
  		theme(legend.position="bottom")
  	return(p)	

}

make_per_tissue_tau_barplot_stratefied_by_annotation_model <- function(df, titler, ordered_annotation_models) {
	frac_lb_arr <- c(df$tau - 1.96*df$tau_se)
	frac_ub_arr <- c(df$tau + 1.96*df$tau_se)

	df2 <- data.frame(tissue=df$tissue,tau=df$tau, tau_lb=frac_lb_arr, tau_ub=frac_ub_arr, method=df$annotation_model)

	df2$method = factor(df2$method, levels=ordered_annotation_models)

	p<-ggplot(data=df2, aes(x=tissue, y=tau, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=tau_lb, ymax=tau_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=titler, x="", y="TGFM-SLDSC Tau", fill="") +
  		theme(legend.position="bottom")

  	return(p)	

}


make_per_tissue_tau_barplot_stratefied_by_sparse_models <- function(df, titler, sparse_model_name, reg_weights) {
	frac_lb_arr <- c(df$tau - 1.96*df$tau_se)
	frac_ub_arr <- c(df$tau + 1.96*df$tau_se)



	df2 <- data.frame(tissue=df$tissue,tau=df$tau, tau_lb=frac_lb_arr, tau_ub=frac_ub_arr, method=df$model_type)

	ordered_methods <- c("tgfm-sldsc")
	for (reg_iter in 1:length(reg_weights)) {
		reg_weight = reg_weights[reg_iter]
		ordered_methods <- c(ordered_methods, paste0(sparse_model_name, "_", reg_weight))
	}


	df2$method = factor(df2$method, levels=ordered_methods)

	p<-ggplot(data=df2, aes(x=tissue, y=tau, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=tau_lb, ymax=tau_ub), width=.2, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
  		labs(title=titler, x="", y="TGFM-SLDSC Tau", fill="") +
  		theme(legend.position="bottom")+guides(fill=guide_legend(nrow=2,byrow=TRUE))

  	return(p)	

}

make_coef_heatmap <- function(df) {

	df$tau[df$tau < 0.0] = 0.0



	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$trait <- recode(df$trait, pigment_HAIR="Hair pigment", body_BALDING1="Balding1",bmd_HEEL_TSCOREz="Heel t-score", biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")
	p <- ggplot(df, aes(x = tissue, y = trait, fill = tau)) +
  		geom_tile() +
  		theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="tau",x="Tissue", y="GWAS trait",title="")
  	return(p)

}


make_coef_z_score_heatmap <- function(df) {
	df$tau_z = df$tau/df$tau_se

	#df$tau_z[(df$tau_z < 2.5) & (df$tau_z > -2.5)] = 0.0
	df$tau_z[df$tau_z < 2.4] = 0.0

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$trait <- recode(df$trait, pigment_HAIR="Hair pigment", body_BALDING1="Balding1",bmd_HEEL_TSCOREz="Heel t-score", biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")
	p <- ggplot(df, aes(x = tissue, y = trait, fill = tau_z)) +
  		geom_tile() +
  		theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low="grey",high="blue") +
  		labs(fill="tau z-score",x="Tissue", y="GWAS trait",title="")
  	return(p)

}

load_in_per_tissue_sparse_tau_df <- function(trait_names, tissue_version, gene_version, annotation_version, regularization_weight, sparse_model_name, tgfm_sldsc_results_dir) {
	tissue_vec <- c()
	trait_vec <- c()
	tau_vec <- c()
	tau_se_vec <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		sparse_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", annotation_version, "_", gene_version, "_", tissue_version, "_pmces_gene_adj_ld_scores_organized_", regularization_weight, "_sparse_", sparse_model_name, "_res.txt")
		t_df <- read.table(sparse_file, header=TRUE, sep="\t")
		first_eqtl_row_num = which(as.character(t_df$Annotation)=="Adipose_Subcutaneous")
		eqtl_df <- t_df[first_eqtl_row_num:dim(t_df)[1],]
		n_tiss = dim(eqtl_df)[1]
		tissue_vec <- c(tissue_vec, as.character(eqtl_df$Annotation))
		tau_vec <- c(tau_vec, eqtl_df$tau)
		tau_se_vec <- c(tau_se_vec, eqtl_df$tau_se)
		trait_vec <- c(trait_vec, rep(trait_name, n_tiss))
	}
	df <- data.frame(trait=trait_vec, tissue=tissue_vec, tau=tau_vec, tau_se=tau_se_vec)
	return(df)
}


independent_trait_names_file = args[1]
tgfm_sldsc_results_dir = args[2]
preprocessed_tgfm_sldsc_data_dir = args[3]
gtex_tissue_colors_file = args[4]
visualize_tgfm_sldsc_dir = args[5]


# Extract trait names
trait_df <- read.table(independent_trait_names_file, header=TRUE, sep="\t")
trait_names <- as.character(trait_df$study_name)


# Models to iterate over
## 1. all_tissues vs non_sex_tissues
## 2. component_gene vs cis_heritable_gene
## 3. Annotation models: 'genotype_intercept', 'baseline_no_qtl' 'baselineLD_no_qtl' 


############################################
# Total expression mediated h2 across traits
#############################################
if (FALSE) {
# Load in data
tissue_versions <- c("all_tissues", "non_sex_tissues")
gene_versions <- c("component_gene", "cis_heritable_gene")
annotation_versions <- c("genotype_intercept", "baseline_no_qtl", "LDanno_only", "baselineLD_no_qtl")
all_trait_names <- c(trait_names, "pigment_HAIR", "body_BALDING1", "bmd_HEEL_TSCOREz")
total_expression_mediated_h2_df <- load_in_total_expression_mediated_h2_df(all_trait_names, tissue_versions, gene_versions, annotation_versions, tgfm_sldsc_results_dir)

# Make barplot showing total expression-mediated heritabilities across traits for annotation models (seperate for each other parameters)
for (tissue_version_iter in 1:length(tissue_versions)) {
	for (gene_version_iter in 1:length(gene_versions)) {
		# Get data for specific gene version for specific tissue version
		tissue_version <- tissue_versions[tissue_version_iter]
		gene_version <- gene_versions[gene_version_iter]

		subset_df <- total_expression_mediated_h2_df[as.character(total_expression_mediated_h2_df$tissue_version)==tissue_version,]
		subset_df <- subset_df[as.character(subset_df$gene_version)==gene_version,]

		plot_title=paste0(as.character(tissue_version), " / ", gene_version)
		p <- make_fraction_mediated_h2_barplot_stratefied_by_annotation_models(subset_df, plot_title, annotation_versions)
		output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_total_expression_mediated_h2_annotation_stratefied_barplot_", tissue_version, "_", gene_version,".pdf")
		ggsave(p, file=output_file, width=7.2, height=4.5, units="in")
	}
}


# Make barplot showing total expression-mediated heritabilities across traits for different gene and tissue models (seperate for each annotation model)
annotation_version <- "baseline_no_qtl"

for (annotation_version_iter in 1:length(annotation_versions)) {
	annotation_version <- annotation_versions[annotation_version_iter]
	subset_df <-total_expression_mediated_h2_df[as.character(total_expression_mediated_h2_df$annotation_model)==annotation_version,]

	plot_title=paste0(annotation_version)
	p <- make_fraction_mediated_h2_barplot_stratefied_by_gene_and_tissue_models(subset_df, plot_title)
	output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_total_expression_mediated_h2_gene_and_tissue_stratefied_barplot_", annotation_version,".pdf")
	ggsave(p, file=output_file, width=7.2, height=4.5, units="in")
}






############################################
# Per-tissue tau estimates
#############################################
# Load in data
tissue_versions <- c("all_tissues", "non_sex_tissues")
gene_versions <- c("component_gene", "cis_heritable_gene")
annotation_versions <- c("genotype_intercept", "baseline_no_qtl", "baselineLD_no_qtl")
gene_tau_df <- load_in_per_tissue_tau_df(trait_names, tissue_versions, gene_versions, annotation_versions, tgfm_sldsc_results_dir)


subset_df <- gene_tau_df[as.character(gene_tau_df$tissue_version)=="non_sex_tissues",]
subset_df <- subset_df[as.character(subset_df$gene_version)=="component_gene",]
subset_df <- subset_df[as.character(subset_df$annotation_model)=="baseline_no_qtl",]





# Make barplot showing per-tissue tau stratefied by annotation model for (seperate for each other parameters and traits)
for (tissue_version_iter in 1:length(tissue_versions)) {
	for (gene_version_iter in 1:length(gene_versions)) {
		# Get data for specific gene version for specific tissue version
		tissue_version <- tissue_versions[tissue_version_iter]
		gene_version <- gene_versions[gene_version_iter]

		subset_df <- gene_tau_df[as.character(gene_tau_df$tissue_version)==tissue_version,]
		subset_df <- subset_df[as.character(subset_df$gene_version)==gene_version,]

		for (trait_iter in 1:length(trait_names)) {
			trait_name <- trait_names[trait_iter]
			
			subset_trait_df <- subset_df[as.character(subset_df$trait)==trait_name,]

			plot_title=paste0(as.character(trait_name, ":    ", tissue_version), " / ", gene_version)
			p <- make_per_tissue_tau_barplot_stratefied_by_annotation_model(subset_trait_df, plot_title, annotation_versions)
			output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_per_tissue_tau_annotation_stratefied_barplot_", trait_name, "_", tissue_version, "_", gene_version,".pdf")

			ggsave(p, file=output_file, width=7.2, height=4.5, units="in")
		}
	}
}

# Make barplot showing per-tissue tau stratefied by different gene and tissue models for (seperate for each other parameters and traits)
tissue_version <- "non_sex_tissues"
for (annotation_version_iter in 1:length(annotation_versions)) {
	for (trait_iter in 1:length(trait_names)) {
		annotation_version <- annotation_versions[annotation_version_iter]
		subset_df <-gene_tau_df[as.character(gene_tau_df$annotation_model)==annotation_version,]

		trait_name <- trait_names[trait_iter]
		subset_df <- subset_df[as.character(subset_df$trait)==trait_name,]
		subset_trait_df <- subset_df[as.character(subset_df$tissue_version)==tissue_version,]

		plot_title=paste0(trait_name, " ", annotation_version)
		p <- make_per_tissue_tau_barplot_stratefied_by_gene_models(subset_trait_df, plot_title)
		output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_per_tissue_tau_gene_stratefied_barplot_", trait_name, "_", annotation_version,".pdf")
		ggsave(p, file=output_file, width=7.2, height=4.5, units="in")
	}
}





models <- c("ard_all_coefficients_mv_update", "ard_eqtl_coefficients_mv_update")
regularization_weights <- c("0.1","0.5", "1.0")
annotation_version <- "baseline_no_qtl"
gene_version <- "component_gene"
tissue_version <- "non_sex_tissues"


for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]

	for (model_iter in 1:length(models)) {
		sparse_model_name <- models[model_iter]
		tau_df <- load_in_per_tissue_tau_df_for_sparse_and_non_sparse_models(trait_name, tissue_version, gene_version, annotation_version, sparse_model_name, regularization_weights, tgfm_sldsc_results_dir)

		plot_title=paste0(trait_name, " ", sparse_model_name)
		p <- make_per_tissue_tau_barplot_stratefied_by_sparse_models(tau_df, plot_title, sparse_model_name, regularization_weights)
		output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_per_tissue_tau_sparse_stratefied_barplot_", trait_name, "_", sparse_model_name, "_", annotation_version, "_", gene_version, "_", tissue_version,".pdf")
		ggsave(p, file=output_file, width=7.2, height=4.5, units="in")		

	}

}


}


############################################
# Per-tissue tau estimates for single model 

#############################################
# Load in data for sparse modell
tissue_version <- "non_sex_tissues"
gene_version <- "component_gene"
annotation_version <- "baseline_no_qtl"
reg_param <- "0.5"
sparse_model_name <- "ard_eqtl_coefficients_mv_update"
sparse_gene_tau_df <- load_in_per_tissue_sparse_tau_df(trait_names, tissue_version, gene_version, annotation_version, reg_param, sparse_model_name, tgfm_sldsc_results_dir)

tau_heatmap <- make_coef_heatmap(sparse_gene_tau_df) + theme(legend.position="right")
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_sparse_coef_heatmap_", tissue_version, "_", gene_version, "_", annotation_version, "_", reg_param,".pdf")
ggsave(tau_heatmap, file=output_file, width=8.0, height=7.5, units="in")	

tissue_version <- "all_tissues"
gene_version <- "component_gene"
annotation_version <- "baseline_no_qtl"
reg_param <- "0.5"
sparse_model_name <- "ard_eqtl_coefficients_mv_update"
sparse_gene_tau_df <- load_in_per_tissue_sparse_tau_df(trait_names, tissue_version, gene_version, annotation_version, reg_param, sparse_model_name, tgfm_sldsc_results_dir)

tau_heatmap <- make_coef_heatmap(sparse_gene_tau_df) + theme(legend.position="right")
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_sparse_coef_heatmap_", tissue_version, "_", gene_version, "_", annotation_version, "_", reg_param,".pdf")
ggsave(tau_heatmap, file=output_file, width=8.0, height=7.5, units="in")	


# Load in data for sparse modell
tissue_version <- "non_sex_tissues"
gene_version <- "component_gene"
annotation_version <- "baselineLD_no_qtl"
reg_param <- "0.5"
sparse_model_name <- "ard_eqtl_coefficients_mv_update"
sparse_gene_tau_df <- load_in_per_tissue_sparse_tau_df(trait_names, tissue_version, gene_version, annotation_version, reg_param, sparse_model_name, tgfm_sldsc_results_dir)

tau_heatmap <- make_coef_heatmap(sparse_gene_tau_df) + theme(legend.position="right")
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_sparse_coef_heatmap_", tissue_version, "_", gene_version, "_", annotation_version, "_", reg_param,".pdf")
ggsave(tau_heatmap, file=output_file, width=8.0, height=7.5, units="in")	

tissue_version <- "all_tissues"
gene_version <- "component_gene"
annotation_version <- "baselineLD_no_qtl"
reg_param <- "0.5"
sparse_model_name <- "ard_eqtl_coefficients_mv_update"
sparse_gene_tau_df <- load_in_per_tissue_sparse_tau_df(trait_names, tissue_version, gene_version, annotation_version, reg_param, sparse_model_name, tgfm_sldsc_results_dir)

tau_heatmap <- make_coef_heatmap(sparse_gene_tau_df) + theme(legend.position="right")
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_sparse_coef_heatmap_", tissue_version, "_", gene_version, "_", annotation_version, "_", reg_param,".pdf")
ggsave(tau_heatmap, file=output_file, width=8.0, height=7.5, units="in")	






#############################################
# Load in data for non-sparse modell
tissue_versions <- c("all_tissues", "non_sex_tissues")
gene_versions <- c("component_gene", "cis_heritable_gene")
annotation_versions <- c("genotype_intercept", "baseline_no_qtl", "baselineLD_no_qtl")
gene_tau_df <- load_in_per_tissue_tau_df(trait_names, tissue_versions, gene_versions, annotation_versions, tgfm_sldsc_results_dir)

# Filter to single model
tissue_version <- "non_sex_tissues"
gene_version <- "component_gene"
annotation_version <- "baseline_no_qtl"
gene_tau_df2 <- gene_tau_df[as.character(gene_tau_df$tissue_version)==tissue_version,]
gene_tau_df2 <- gene_tau_df2[as.character(gene_tau_df2$gene_version)==gene_version,]
gene_tau_df2 <- gene_tau_df2[as.character(gene_tau_df2$annotation_model)==annotation_version,]


# Make heatmap showing z-scores of multivariable approach across traits
z_score_heatmap <- make_coef_z_score_heatmap(gene_tau_df2)
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_coef_z_score_heatmap_", tissue_version, "_", gene_version, "_", annotation_version,".pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=7.5, units="in")		


# Filter to single model
tissue_version <- "all_tissues"
gene_version <- "component_gene"
annotation_version <- "baseline_no_qtl"
gene_tau_df2 <- gene_tau_df[as.character(gene_tau_df$tissue_version)==tissue_version,]
gene_tau_df2 <- gene_tau_df2[as.character(gene_tau_df2$gene_version)==gene_version,]
gene_tau_df2 <- gene_tau_df2[as.character(gene_tau_df2$annotation_model)==annotation_version,]


# Make heatmap showing z-scores of multivariable approach across traits
z_score_heatmap <- make_coef_z_score_heatmap(gene_tau_df2)
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_coef_z_score_heatmap_", tissue_version, "_", gene_version, "_", annotation_version,".pdf")
ggsave(z_score_heatmap, file=output_file, width=7.2, height=7.5, units="in")		

















































if (FALSE) {
	heatmap <- make_ld_score_annotation_correlation_heatmap(preprocessed_tgfm_sldsc_data_dir)
	output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_annotation_correlation_heatmap.pdf")
	ggsave(heatmap, file=output_file, width=9.2, height=10.5, units="in")
}





if (FALSE) {
for (trait_iter in 1:length(trait_names)) {
trait_name <- trait_names[trait_iter]
# Extract sparse model df
models <- c("ard_all_coefficients_mv_update")
regularization_weights <- c("0.1", "1.0", "2.0", "5.0")
tau_df_all_coef <- load_in_tgfm_sldsc_coefficients(tgfm_sldsc_results_dir, trait_name, models, regularization_weights)
models <- c("ard_eqtl_coefficients_mv_update")
regularization_weights <- c("0.1", "1.0", "2.0", "5.0")
tau_df_eqtl_coef <- load_in_tgfm_sldsc_coefficients(tgfm_sldsc_results_dir, trait_name, models, regularization_weights)


# For this trait make boxplot showing sparse effect size
model_name <- "ard_all_coefficients_mv_update"
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_",model_name, "_per_tissue_tau_for_range_of_regularizations_boxplot_", trait_name, ".pdf")
barplot <- per_tissue_barplot_for_range_of_regularizations(tau_df_all_coef, model_name, regularization_weights, trait_name)
ggsave(barplot, file=output_file, width=7.2, height=4.5, units="in")

model_name <- "ard_eqtl_coefficients_mv_update"
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_",model_name, "_per_tissue_tau_for_range_of_regularizations_boxplot_", trait_name, ".pdf")
barplot <- per_tissue_barplot_for_range_of_regularizations(tau_df_eqtl_coef, model_name, regularization_weights, trait_name)
ggsave(barplot, file=output_file, width=7.2, height=4.5, units="in")
}


for (trait_iter in 1:length(trait_names)) {
trait_name <- trait_names[trait_iter]
# Extract sparse model df
models <- c("ard_all_coefficients")
regularization_weights <- c("1.0", "2.0", "5.0", "10.0")
tau_df_all_coef <- load_in_tgfm_sldsc_coefficients(tgfm_sldsc_results_dir, trait_name, models, regularization_weights)
models <- c("ard_eqtl_coefficients")
regularization_weights <- c("1.0", "2.0", "5.0", "10.0")
tau_df_eqtl_coef <- load_in_tgfm_sldsc_coefficients(tgfm_sldsc_results_dir, trait_name, models, regularization_weights)


# For this trait make boxplot showing sparse effect size
model_name <- "ard_all_coefficients"
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_",model_name, "_per_tissue_tau_for_range_of_regularizations_boxplot_", trait_name, ".pdf")
barplot <- per_tissue_barplot_for_range_of_regularizations(tau_df_all_coef, model_name, regularization_weights, trait_name)
ggsave(barplot, file=output_file, width=7.2, height=4.5, units="in")

model_name <- "ard_eqtl_coefficients"
output_file <- paste0(visualize_tgfm_sldsc_dir, "tgfm_sldsc_",model_name, "_per_tissue_tau_for_range_of_regularizations_boxplot_", trait_name, ".pdf")
barplot <- per_tissue_barplot_for_range_of_regularizations(tau_df_eqtl_coef, model_name, regularization_weights, trait_name)
ggsave(barplot, file=output_file, width=7.2, height=4.5, units="in")
}

}
