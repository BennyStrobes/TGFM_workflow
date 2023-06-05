args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(ggbeeswarm)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

extract_data_frame_containing_fraction_of_mediated_components_in_each_tissue_per_trait <- function(trait_names, tgfm_results_dir) {
	trait_names_arr <- c()
	tissue_names_arr <- c()
	fraction_mediated_arr <- c()
	expr_mediated_fraction_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		tiss_summary_file <- paste0(tgfm_results_dir, trait_name, "_tgfm_component_tissue_prob_summary.txt")
		tiss_data <- read.table(tiss_summary_file, header=TRUE, sep="\t")
		tiss_names <- colnames(tiss_data)[4:length(colnames(tiss_data))]
		tiss_data = tiss_data[,4:(dim(tiss_data)[2])]
		average_med_prob_per_tissue = as.numeric(colSums(tiss_data))/dim(tiss_data)[1]
		tissue_names_arr <- c(tissue_names_arr, tiss_names)
		fraction_mediated_arr <- c(fraction_mediated_arr, average_med_prob_per_tissue)
		expr_mediated_fraction_arr <- c(expr_mediated_fraction_arr, average_med_prob_per_tissue/sum(average_med_prob_per_tissue))
		trait_names_arr <- c(trait_names_arr, rep(trait_name, length(average_med_prob_per_tissue)))
	}	
	df = data.frame(trait=trait_names_arr, tissue=tissue_names_arr, average_mediated_probability=fraction_mediated_arr, fraction_expr_mediated_probability=expr_mediated_fraction_arr)
	return(df)
}

extract_data_frame_containing_fraction_of_mediated_components_per_trait <- function(trait_names, tgfm_results_dir) {
	avg_med_prob_arr <- c()
	se_avg_med_prob_arr <- c()
	n_components_arr <- c()
	trait_names_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		tiss_summary_file <- paste0(tgfm_results_dir, trait_name, "_tgfm_predicted_causal_effect_size.txt")
		aa <- read.table(tiss_summary_file, header=TRUE)
		probs = aa$gene_mediated_probability

		avg_med_prob_arr <- c(avg_med_prob_arr, mean(probs))
		se_avg_med_prob_arr <- c(se_avg_med_prob_arr, std_mean(probs))
		n_components_arr <- c(n_components_arr, length(probs))
		trait_names_arr <- c(trait_names_arr, trait_name)
	}

	df <- data.frame(average_mediated_probability=avg_med_prob_arr, average_mediated_probability_se=se_avg_med_prob_arr, trait=trait_names_arr, n_components=n_components_arr)
	return(df)
}

extract_expression_mediated_probabilities_across_components <- function(trait_name, tgfm_results_dir) {
	tiss_summary_file <- paste0(tgfm_results_dir, trait_name, "_tgfm_predicted_causal_effect_size.txt")
	aa <- read.table(tiss_summary_file, header=TRUE)
	probs = aa$gene_mediated_probability
	return(probs)
}


std_mean <- function(x) sd(x)/sqrt(length(x))

make_fraction_of_genetic_elements_from_gene_expression_se_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir) {

	trait_names_vec <- c()
	fraction_mediated_vec <- c()
	fraction_mediated_se_vec <- c()
	pip_threshold_vec <- c()

	pip_thresholds <- c("0.25", "0.5", "0.75")


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]

		for (pip_iter in 1:length(pip_thresholds)) {
			pip_threshold <- pip_thresholds[pip_iter]

			summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_threshold, ".txt")
			tmp <- read.table(summary_file, header=TRUE)


			variant_count = tmp$count[1]
			gene_count = tmp$count[2]
			total_count = gene_count + variant_count

			mean_mediated_prob <- gene_count/(total_count)
		
			se_mediated_prob = sqrt((mean_mediated_prob*(1.0-mean_mediated_prob))/(total_count))
			trait_names_vec <- c(trait_names_vec, trait_name_readable)
			fraction_mediated_vec <- c(fraction_mediated_vec, mean_mediated_prob)
			fraction_mediated_se_vec <- c(fraction_mediated_se_vec, se_mediated_prob)
			pip_threshold_vec <- c(pip_threshold_vec, pip_threshold)
		}
	}

	df <- data.frame(trait=trait_names_vec, average_mediated_probability=fraction_mediated_vec, average_mediated_probability_se=fraction_mediated_se_vec, pip=factor(pip_threshold_vec,levels=pip_thresholds))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$pip) == "0.25",]

	indices = order(tmp_df$average_mediated_probability)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=average_mediated_probability, fill=pip)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		scale_fill_manual(values=c("plum2", "orchid2", "magenta2")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    		labs(y="Fraction of expression-mediated\ngenetic elements", x="") +
    		geom_errorbar( aes(ymin=average_mediated_probability-(1.96*average_mediated_probability_se), ymax=average_mediated_probability+(1.96*average_mediated_probability_se)), width=0.2, colour="grey50", position=position_dodge(.9)) +
    		figure_theme()
    return(p)
}

make_mediated_prob_se_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir) {

	trait_names_vec <- c()
	fraction_mediated_vec <- c()
	fraction_mediated_se_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]

		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_component_level_summary.txt")
		
		tmp <- read.table(summary_file, header=TRUE)
		mediated_probs <- 1.0 - tmp$non_mediated_probability
		
		mean_mediated_prob = mean(mediated_probs)
		se_mediated_prob = sqrt((mean_mediated_prob*(1.0-mean_mediated_prob))/(length(mediated_probs)/100.0))
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		fraction_mediated_vec <- c(fraction_mediated_vec, mean_mediated_prob)
		fraction_mediated_se_vec <- c(fraction_mediated_se_vec, se_mediated_prob)
	}

	df <- data.frame(trait=trait_names_vec, average_mediated_probability=fraction_mediated_vec, average_mediated_probability_se=fraction_mediated_se_vec)


	#df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	indices = order(df$average_mediated_probability)
	df$trait = factor(df$trait, levels=as.character(df$trait)[indices])

	p <- ggplot(df) +
    		geom_bar( aes(x=trait, y=average_mediated_probability), stat="identity", fill="skyblue", alpha=0.7) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    		labs(y="Average expression-mediated\nprobability", x="") +
    		geom_errorbar( aes(x=trait, ymin=average_mediated_probability-(1.96*average_mediated_probability_se), ymax=average_mediated_probability+(1.96*average_mediated_probability_se)), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    		figure_theme()
    return(p)
}

 extract_anno_size <- function(processed_tgfm_sldsc_data_dir) {
 	for (chrom_num in 1:22) {
 		chrom_filer = paste0(processed_tgfm_sldsc_data_dir, "baselineLD_no_qtl_pmces_gene_adj_ld_scores.", chrom_num, ".l2.M_5_50")
 		tmp_data = read.table(chrom_filer,header=FALSE)
 		if (chrom_num == 1) {
 			counter = as.numeric(tmp_data[1,])
 		} else {
 			counter = counter + as.numeric(tmp_data[1,])
 		}
 	}
 	return(counter)
}

extract_data_frame_containing_fraction_of_mediated_h2_in_each_tissue_per_trait <- function(full_trait_names, tgfm_sldsc_results_dir, methods, anno_size, nonneg=TRUE) {
	trait_arr <- c()
	tissue_arr <- c()
	h2_med_arr <- c()
	stand_h2_med_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]
		method <- methods[1]
		per_ele_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baselineLD_no_qtl_pmces_gene_adj_ld_scores_organized_", method, "_res.txt")
		tmp_df <- read.table(per_ele_h2_file, header=TRUE)

		partition_h2 = tmp_df$tau*anno_size
		geno_h2 = sum(partition_h2[1:93])
		tiss_names = as.character(tmp_df$Annotation[94:length(partition_h2)])
		per_tissue_h2 = partition_h2[94:length(partition_h2)]
		if (nonneg==TRUE) {
			per_tissue_h2[per_tissue_h2 < 0.0] = 0.0
		}


		trait_arr <- c(trait_arr, rep(trait_name, length(tiss_names)))
		tissue_arr <- c(tissue_arr, tiss_names)
		h2_med_arr <- c(h2_med_arr, per_tissue_h2)
		stand_h2_med_arr <- c(stand_h2_med_arr, per_tissue_h2/sum(per_tissue_h2))
	}

	df = data.frame(trait=trait_arr, tissue=tissue_arr, average_mediated_probability=h2_med_arr, fraction_expr_mediated_probability=stand_h2_med_arr)

	return(df)
}

extract_fraction_of_h2_mediated_by_gene_expression_for_sparse_models <- function(full_trait_names, tgfm_sldsc_results_dir, methods, anno_size, nonneg=FALSE) {
	trait_arr <- c()
	method_arr <- c()
	med_h2_arr <- c()
	med_h2_se_arr <- c()
	med_h2_lb_arr <- c()
	med_h2_ub_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]
		#hr_trait_name <- hr_trait_names[trait_iter]

		for (method_iter in 1:length(methods)) {
			method <- methods[method_iter]
			per_ele_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baselineLD_no_qtl_pmces_gene_adj_ld_scores_organized_", method, "_res.txt")
			tmp_df <- read.table(per_ele_h2_file, header=TRUE)

			partition_h2 = tmp_df$tau*anno_size
			geno_h2 = sum(partition_h2[1:93])

			per_tissue_h2 = partition_h2[94:length(partition_h2)]
			if (nonneg==TRUE) {
				per_tissue_h2[per_tissue_h2 < 0.0] = 0.0
			}
			#print(per_tissue_h2)

			expr_h2 = sum(per_tissue_h2)
	
			med_h2 = expr_h2/(expr_h2 + geno_h2)

			trait_arr <- c(trait_arr, trait_name)
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


extract_fraction_of_h2_mediated_by_gene_expression_for_several_methods <- function(full_trait_names, tgfm_sldsc_results_dir, methods) {
	trait_arr <- c()
	method_arr <- c()
	med_h2_arr <- c()
	med_h2_se_arr <- c()
	med_h2_lb_arr <- c()
	med_h2_ub_arr <- c()
	for (trait_iter in 1:length(full_trait_names)) {
		trait_name <- full_trait_names[trait_iter]

		for (method_iter in 1:length(methods)) {
			method <- methods[method_iter]
			expr_med_h2_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_", method, "_gene_adj_ld_scores_h2_5_50_med.txt")
			tmp_df <- read.table(expr_med_h2_file, header=TRUE)

			trait_arr <- c(trait_arr, trait_name)
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

make_scatterplot_comparing_expression_mediated_h2_across_traits <- function(expr_med_frac_df, expr_med_frac_sparse_model_df, x_axis_label, y_axis_label) {
	df <- data.frame(expr_mediated_v1=expr_med_frac_df$fraction_h2, expr_mediated_v2=expr_med_frac_sparse_model_df$fraction_h2)
	corry=cor(expr_med_frac_df$fraction_h2, expr_med_frac_sparse_model_df$fraction_h2)
	p <- ggplot(df, aes(x=expr_mediated_v1, y=expr_mediated_v2)) + geom_point() +
	figure_theme() + 
	geom_abline()+
	labs(x=x_axis_label, y=y_axis_label,title=paste0('correlation: ', corry))
	return(p)
}
make_scatterplot_comparing_expression_mediated_h2_and_average_component_mediated_probability_across_traits <- function(expr_med_frac_df, expr_med_frac_sparse_model_df, x_axis_label, y_axis_label) {
	df <- data.frame(expr_mediated_v1=expr_med_frac_df$fraction_h2, expr_mediated_v2=expr_med_frac_sparse_model_df$average_mediated_probability)
	corry=cor(expr_med_frac_df$fraction_h2, expr_med_frac_sparse_model_df$average_mediated_probability)
	p <- ggplot(df, aes(x=expr_mediated_v1, y=expr_mediated_v2)) + geom_point() +
	figure_theme() + 
	geom_abline()+
	labs(x=x_axis_label, y=y_axis_label,title=paste0('correlation: ', corry))
	return(p)
}

plot_distribution_of_tgfm_expression_mediated_probabilities <- function(expression_med_prob,titler) {
	df <- data.frame(tgfm_expression_mediated_probability=expression_med_prob)
	p <- ggplot(df, aes(x=tgfm_expression_mediated_probability))+
  		geom_histogram(color="darkblue", fill="lightblue") +
  		figure_theme() +
  		labs(x="TGFM expression-mediated probability", title=titler)
}

stacked_barplot_for_single_trait_breaking_down_expression_effects_per_tissue <- function(trait_tissue_fraction_component_mediated_df, gtex_colors_df, trait_name,x_axis_label) {
	df <- data.frame(tissue=as.character(trait_tissue_fraction_component_mediated_df$tissue), fraction=trait_tissue_fraction_component_mediated_df$fraction_expr_mediated_probability)

	tissue_arr = as.character(df$tissue)
	fraction_arr = df$fraction
	# Summarize by top 3 tissues
	ordered_tissues <- tissue_arr[order(-fraction_arr)]

	tissue_arr2 <- c()
	fraction_arr2 <- c()
	trait_arr <- c()
	colors <- c()
	num_tissues <- 4
	for (indexer in 1:num_tissues) {
		if (ordered_tissues[indexer] == "Non-mediated") {
			colors <- c(colors, "#999999")
		} else {
			tt <- ordered_tissues[indexer]
			if (startsWith(tt, "Brain")) {
				colors <- c(colors, "#EEEE00")
			} 
			else if (startsWith(tt, "Artery_Heart")) {
				colors <- c(colors, "#FF5555")
			}
			else {
				index = which(gtex_colors_df$tissue_site_detail_id==tt)
				colors <- c(colors, paste0('#', as.character(gtex_colors_df$tissue_color_hex)[index]))
			}
		}

		tissue_arr2 <- c(tissue_arr2, paste0(ordered_tissues[indexer], "-mediated"))
		fraction_arr2 <- c(fraction_arr2, as.numeric(df$fraction[df$tissue == ordered_tissues[indexer]]))
		trait_arr <- c(trait_arr, trait_name)
	}

	tissue_arr2 <- c(tissue_arr2, "Other-expression-mediated")
	fraction_arr2 <- c(fraction_arr2, 1.0 - sum(fraction_arr2))
	trait_arr <- c(trait_arr, trait_name)
	colors <- c(colors, "darkblue")

	df2 <- data.frame(tissue=factor(tissue_arr2, levels=rev(c(tissue_arr2))), fraction=fraction_arr2, trait=trait_arr, color=colors)


	p <- ggplot(df2, aes(fill=tissue, y=trait, x=fraction)) + 
    	geom_bar(position="fill", stat="identity") +
    	figure_theme() +
    	labs(x=x_axis_label, title=trait_name, y="", fill="") +
    	scale_fill_manual(values=rev(as.character(df2$color))) +
    	theme(legend.position="bottom") +
    	theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank())  +
    	guides(fill = guide_legend(reverse=TRUE, nrow=3, byrow=TRUE)) +
    	theme(legend.key.size = unit(.5, 'cm'), legend.title = element_text(size=11),legend.text = element_text(size=11)) +
    	theme(plot.title = element_text(hjust = 0.5))
   	return(p)

}

make_bar_plot_showing_expected_number_of_causal_tissue_categories_for_single_trait <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir) {
	tissue_vec <- c()
	count_vec <- c()
	pip_thresh_vec <- c()


	pip_threshs <- c(.3,.5,.7)
	trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_tissue_category_pip_summary.txt")
	trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")

	ordered_tissues <- sort(unique(as.character(trait_df$tissue_category_name)))

	tmper <- c()
	for (pip_iter in 1:length(pip_threshs)) {
		pip_thresh <- pip_threshs[pip_iter]
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		#print(head(tmp_df))
		for (tissue_iter in 1:length(ordered_tissues)) {
			tissue_name <- ordered_tissues[tissue_iter]
			pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_category_name) ==tissue_name])
			tissue_vec <- c(tissue_vec, tissue_name)
			count_vec <- c(count_vec, pip_sum)
			pip_thresh_vec <- c(pip_thresh_vec, pip_thresh)
			if (pip_thresh == pip_threshs[1]) {
				tmper <- c(tmper, pip_sum)
			}
		}


	}

	ord <- order(tmper)
	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues[ord]), pip_thresh=factor(pip_thresh_vec))
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	#df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	ordered_tissues2 <- as.character(df$tissue)[1:length(unique(df$tissue))]
	df$tissue = factor(df$tissue, levels=ordered_tissues2[ord])

	p<-ggplot(df, aes(x=tissue, y=expected_causal_genes, fill=pip_thresh)) +
  		geom_bar(stat="identity",position=position_dodge())+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="bottom") +
  		labs(x="", y="Expected #\ncausal tissues", fill="PIP threshold",title=trait_name_readable)

  	return(p)
}

make_bar_plot_showing_sldsc_tau_of_each_tissue <- function(sldsc_file, trait_name_readable) {
	# Load in input data
	df <- read.table(sldsc_file, header=TRUE)
	# Skip variants
	df <- df[54:(dim(df)[1]),]

	df$tissue = df$Annotation
	ordered_tissues <- as.character(df$tissue)
	df$prior = df$tau
	priors <- as.numeric(df$prior)

	ord <- order(priors)
	df$tissue = factor(df$tissue, levels=ordered_tissues[ord])
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	ordered_tissues2 <- as.character(df$tissue)[1:length(unique(df$tissue))]
	df$tissue = factor(df$tissue, levels=ordered_tissues2[ord])
	
	p<-ggplot(df, aes(x=tissue, y=priors)) +
  		geom_bar(stat="identity",position=position_dodge())+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="bottom") +
  		geom_errorbar( aes(x=tissue, ymin=priors-(1.96*tau_se), ymax=priors+(1.96*tau_se)), width=0.2, colour="blue", alpha=0.6, size=.8) +
  		labs(x="", y="TGLR Tau", title=trait_name_readable)
  	return(p)
}

make_violin_plot_showing_distribution_of_bootstrapped_taus_of_each_tissue <- function(iterative_prior_file, method_version, trait_name_readable, threshold=1e-4) {
	# Load in input data
	df <- read.table(iterative_prior_file, header=TRUE)
	# Skip variants
	nrow = dim(df)[1]

	df <- df[(nrow-36):(dim(df)[1]),]
	tmp_tissue <- c()
	for (tiss_iter in 1:length(df$annotation_name)) {
		tiss_name <- as.character(df$annotation_name[tiss_iter])
		new_tiss_name = strsplit(tiss_name,"_0")[[1]]
		tmp_tissue <- c(tmp_tissue, new_tiss_name)
	}


	df$tissue = as.character(tmp_tissue)

	tissue_vec <- c()
	prior_prob_vec <- c()
	siggy <- c()


	for (tissue_iter in 1:length(df$tissue)) {
		tissue_name <- as.character(df$tissue[tissue_iter])
		prob_string = as.character(df$bootstrapped_taus[tissue_iter])
		prob_vec = as.numeric(strsplit(prob_string,";")[[1]])

		pvalue = sum(prob_vec < threshold)/length(prob_vec)
		if (pvalue < (.1/37)) {
			sigger = 'significant'
		} else {
			sigger = 'non_significant'
		}

		prior_prob_vec <- c(prior_prob_vec, prob_vec)
		tissue_vec <- c(tissue_vec, rep(tissue_name, length(prob_vec)))
		siggy <- c(siggy, rep(sigger, length(prob_vec)))
	}


	df2 <- data.frame(tissue=factor(tissue_vec, levels=as.character(df$tissue)), probability=prior_prob_vec, significance=factor(siggy))
	df2$tissue = str_replace_all(as.character(df2$tissue), "-", "_")
	df2$tissue <- recode(df2$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	df2$tissue = factor(df2$tissue)

    p <- ggplot(df2, aes(tissue, probability)) +
  		geom_boxplot(colour = "grey50", outlier.shape = NA) +
  		geom_jitter(aes(color=siggy),size=.001) +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		labs(x="", y="Per gene-tissue h2", title=paste0(trait_name_readable, " / ", method_version)) +
  		figure_theme() +
  		theme(legend.position="None") +
  		scale_colour_manual(values = c("grey50", "darkorchid3"))


  	return(p)


}


make_violin_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue <- function(iterative_prior_file, method_version, trait_name_readable, threshold=1e-8) {
	# Load in input data
	df <- read.table(iterative_prior_file, header=TRUE)
	# Skip variants
	nrow = dim(df)[1]

	df <- df[(nrow-36):(dim(df)[1]),]
	df$tissue = df$element_name

	tissue_vec <- c()
	prior_prob_vec <- c()
	siggy <- c()


	for (tissue_iter in 1:length(df$tissue)) {
		tissue_name <- as.character(df$tissue[tissue_iter])
		prob_string = as.character(df$prior_distribution[tissue_iter])
		prob_vec = as.numeric(strsplit(prob_string,";")[[1]])

		pvalue = sum(prob_vec < threshold)/length(prob_vec)
		if (pvalue < (.1/37)) {
			sigger = 'significant'
		} else {
			sigger = 'non_significant'
		}

		prior_prob_vec <- c(prior_prob_vec, prob_vec)
		tissue_vec <- c(tissue_vec, rep(tissue_name, length(prob_vec)))
		siggy <- c(siggy, rep(sigger, length(prob_vec)))
	}


	df2 <- data.frame(tissue=factor(tissue_vec, levels=as.character(df$tissue)), probability=prior_prob_vec, significance=factor(siggy))
	df2$tissue = str_replace_all(as.character(df2$tissue), "-", "_")
	df2$tissue <- recode(df2$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	df2$tissue = factor(df2$tissue)

    p <- ggplot(df2, aes(tissue, probability)) +
  		geom_boxplot(colour = "grey50", outlier.shape = NA) +
  		geom_jitter(aes(color=siggy),size=.001) +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		labs(x="", y="Prior Probability", title=paste0(trait_name_readable, " / ", method_version)) +
  		figure_theme() +
  		theme(legend.position="None") +
  		scale_colour_manual(values = c("grey50", "darkorchid3"))


  	return(p)


}

make_bar_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue <- function(iterative_prior_file, method_version, trait_name_readable) {
	# Load in input data
	df <- read.table(iterative_prior_file, header=TRUE)
	# Skip variants
	df <- df[2:(dim(df)[1]),]
	df$tissue = df$element_name

	ordered_tissues <- as.character(df$tissue)
	priors <- as.numeric(df$prior)

	lb_ci = c()
	ub_ci = c()

	nrows = dim(df)[1]
	for (row_iter in 1:nrows) {
		stringer = as.character(df$prior_distribution[row_iter])
		distribution = as.numeric(strsplit(stringer, split=";")[[1]])
		n_samples = length(distribution)
		sides = floor(n_samples*.025)
		sorted_distr = sort(distribution)
		lb_ci <- c(lb_ci, sorted_distr[sides])
		ub_ci <- c(ub_ci, sorted_distr[(n_samples-sides)])
	}
	df$ub_ci = ub_ci
	df$lb_ci = lb_ci

	ord <- order(priors)
	df$tissue = factor(df$tissue, levels=ordered_tissues[ord])
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	ordered_tissues2 <- as.character(df$tissue)[1:length(unique(df$tissue))]
	df$tissue = factor(df$tissue, levels=ordered_tissues2[ord])

	
	p<-ggplot(df, aes(x=tissue, y=priors)) +
  		geom_bar(stat="identity",position=position_dodge())+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="bottom") +
  		labs(x="", y="Prior probability", title=trait_name_readable) +
  		geom_errorbar( aes(x=tissue, ymin=lb_ci, ymax=ub_ci), width=0.2, colour="blue", alpha=0.6, size=.8)

}


make_bar_plot_showing_iterative_prior_probability_of_each_tissue <- function(iterative_prior_file, method_version, trait_name_readable) {
	# Load in input data
	df <- read.table(iterative_prior_file, header=TRUE)
	# Skip variants
	df <- df[2:(dim(df)[1]),]
	df$tissue = df$element_name

	ordered_tissues <- as.character(df$tissue)
	priors <- as.numeric(df$prior)

	ord <- order(priors)
	df$tissue = factor(df$tissue, levels=ordered_tissues[ord])
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	ordered_tissues2 <- as.character(df$tissue)[1:length(unique(df$tissue))]
	df$tissue = factor(df$tissue, levels=ordered_tissues2[ord])
	
	p<-ggplot(df, aes(x=tissue, y=priors)) +
  		geom_bar(stat="identity",position=position_dodge())+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="bottom") +
  		labs(x="", y="Prior probability", title=trait_name_readable)

}

make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir, ordered_tissues) {
	tissue_vec <- c()
	count_vec <- c()
	pip_thresh_vec <- c()


	if (FALSE) {
	if (method_version == "susie_sampler_iterative_variant_gene_tissue") {
		input_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", "susie_pmces_variant_gene", "_tgfm_per_tissue_sldsc_comparison.txt")
	} else {
		input_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_tissue_sldsc_comparison.txt")
	}
	df <- read.table(input_file, header=TRUE, sep="\t")
	ordered_tissues <- as.character(sort(df$tissue))
	}

	pip_threshs <- c(.1, .3, .5)
	trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
	trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")

	#ordered_tissues <- sort(unique(as.character(trait_df$tissue_name)))

	tmper <- c()
	for (pip_iter in 1:length(pip_threshs)) {
		pip_thresh <- pip_threshs[pip_iter]
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		#print(head(tmp_df))
		for (tissue_iter in 1:length(ordered_tissues)) {
			tissue_name <- ordered_tissues[tissue_iter]
			pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
			tissue_vec <- c(tissue_vec, tissue_name)
			count_vec <- c(count_vec, pip_sum)
			pip_thresh_vec <- c(pip_thresh_vec, pip_thresh)
			if (pip_thresh == pip_threshs[1]) {
				tmper <- c(tmper, pip_sum)
			}
		}


	}

	ord <- order(tmper)
	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues[ord]), pip_thresh=factor(pip_thresh_vec))
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	ordered_tissues2 <- as.character(df$tissue)[1:length(unique(df$tissue))]
	df$tissue = factor(df$tissue, levels=ordered_tissues2[ord])

	p<-ggplot(df, aes(x=tissue, y=expected_causal_genes, fill=pip_thresh)) +
  		geom_bar(stat="identity",position=position_dodge())+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="bottom") +
  		labs(x="", y="Expected #\ncausal genes", fill="PIP threshold",title=trait_name_readable)

  	return(p)
}


make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir) {
	# First extract data from each trait
	trait_vec <- c()
	tissue_vec <- c()
	broad_tissue_vec <- c()
	pip_vec <- c()
	# Extract data from each trait
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]

		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")
		trait_df <- trait_df[trait_df$PIP > .4, ]

		pip_vec <- c(pip_vec, trait_df$PIP)
		broad_tissue_vec <- c(broad_tissue_vec, as.character(trait_df$tissue_visualization_category))
		tissue_vec <- c(tissue_vec, as.character(trait_df$tissue_name))
		trait_vec <- c(trait_vec, rep(trait_name_readable, length(trait_df$PIP)))
	}
	# Put in organized df
	df <- data.frame(trait=factor(trait_vec), tissue=factor(tissue_vec), tissue_category=factor(broad_tissue_vec), pip=pip_vec)


	# Make beeswarm plot
	# Beeswarm plot in ggplot2
	p <- ggplot(df, aes(x = trait, y = pip, color = tissue_category)) +
  		geom_beeswarm(cex = 1.7) +
  		scale_color_manual(values=c("#FF6600", "#FF00BB", "#EEEE00", "#660099", "#EEBB77", "#AAFF99", "#552200", "#AAEEFF", "#AABB66", "#99FF00", "#AAAAFF", "grey", "#0000FF")) +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="bottom")
  	return(p)
}


make_scatterplot_comparing_tglr_expression_h2_and_n_tgfm_pip_across_tissues <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir) {
	input_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_tissue_sldsc_comparison.txt")
	df <- read.table(input_file, header=TRUE, sep="\t")
	#print(df)
	p <- ggplot(df, aes(x=n_tgfm_components, y=tglr_mediated_h2)) +
  		geom_point() +
  		figure_theme() + 
  		labs(title=trait_name, x="TGFM expected #\ncausal genes (PIP >=.2)", y="TGLR expected mediated h2")
  	return(p)

}

make_bar_plot_showing_tgfm_tissue_pvalues_for_single_trait <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names) {
	input_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_causal_tissue_pvalues.txt")

	df <- read.table(input_file, header=TRUE)
	df$tissue_name = str_replace_all(as.character(df$tissue_name), "-", "_")
	#tissue_vec = as.character(df$tissue_name)
	df$log_p = -log10(df$pvalue)
	ord <- order(df$log_p)
	df$tissue_name <- recode(df$tissue_name, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	ordered_tissues2 <- as.character(df$tissue_name)[1:length(unique(df$tissue_name))]
	df$tissue_name = factor(df$tissue_name, levels=ordered_tissues2[ord])

	n_tiss = length(unique(ordered_tissues2))

	# Bonferonni p-value correction
	bonf_thresh = .05/n_tiss

	df$significant = df$pvalue <= bonf_thresh


	p<-ggplot(df, aes(x=tissue_name, y=log_p, fill=significant)) +
  		geom_bar(stat="identity",position=position_dodge())+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="none") +
  		labs(x="", y="-log10(p)",title=trait_name_readable) +
  		scale_fill_manual(values = c("grey50", "darkorchid3")) +
  		geom_hline(yintercept=-log10(bonf_thresh), linetype=2, color="grey50")

  	return(p)
}


make_tissue_overlap_jaccard_index_heatmap <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names) {
	jaccard_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_tissue_overlap_jaccard.txt")
	df2 <- read.table(jaccard_file, header=TRUE)

	new_tissue1 = c(as.character(df2$tissue1), as.character(df2$tissue2))
	new_tissue2 = c(as.character(df2$tissue2), as.character(df2$tissue1))
	new_jaccard_index= c(df2$jaccard_index, df2$jaccard_index)

	df <- data.frame(tissue1=new_tissue1, tissue2=new_tissue2, jaccard_index=new_jaccard_index)

	df$tissue1 = str_replace_all(as.character(df$tissue1), "-", "_")
	df$tissue2 = str_replace_all(as.character(df$tissue2), "-", "_")
	df$tissue1 <- recode(df$tissue1, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$tissue2 <- recode(df$tissue2, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")


	p <- ggplot(df, aes(tissue1, tissue2, fill=jaccard_index)) + 
  	geom_tile() +
  	scale_fill_gradient(low="white", high="blue") +
  	figure_theme() +
  	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  	return(p)
}



##################################
# Extract command line arguments
##################################

independent_trait_names_file <- args[1]
tgfm_sldsc_results_dir <- args[2]
tgfm_results_dir <- args[3]
processed_tgfm_sldsc_data_dir <- args[4]
gtex_tissue_colors_file <- args[5]
visualize_tgfm_dir <- args[6]
iterative_tgfm_prior_dir <- args[7]



##########################################################
# Load in data
##########################################################

# Load in gtex tissue colors 
gtex_colors_df <- read.table(gtex_tissue_colors_file, header=TRUE, sep="\t")
gtex_colors_df$tissue_site_detail_id = as.character(gtex_colors_df$tissue_site_detail_id)
gtex_colors_df$tissue_site_detail_id[23] = "Cells_Cultured_fibroblasts"


# Extract trait names
print(independent_trait_names_file)
trait_df <- read.table(independent_trait_names_file, header=TRUE, sep="\t")
trait_names <- as.character(trait_df$study_name)
trait_names_readable <- as.character(trait_df$study_name_readable)
#trait_names <- c("biochemistry_Cholesterol", "biochemistry_VitaminD", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "blood_MEAN_PLATELET_VOL", "blood_MONOCYTE_COUNT", "body_BMIz", "body_WHRadjBMIz", "bp_DIASTOLICadjMEDz", "lung_FEV1FVCzSMOKE")
#trait_names_readable <- c("Cholesterol", "VitaminD", "Reticulocyte_count", "Platelet_vol", "Monocyte_count", "BMI", "WHRadjBMI", "Diastolic_BP", "FEV1FVC")




# Extract tissue names
iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_names[1], "_component_gene_", "susie_pmces_uniform", "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
aa = read.table(iterative_prior_file, header=TRUE,sep="\t")
tissue_names = as.character(aa$element_name[2:length(aa$element_name)])

##################################################
# Violin plot showing bootstrapped prior distributions 
##################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]
	print(trait_name)
	# PMCES approach
	#method_version="susie_pmces_uniform"
	#iterative_prior_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_iterative_variant_gene_prior_v1_pip_level_bootstrapped.txt")
	#iterative_pmces_violin_plot <- make_violin_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)

	# Sampler approach
	method_version="susie_pmces_uniform"
	iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
	iterative_sampler_violin_plot <- make_violin_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)

	# Sampler approach
	#method_version="susie_pmces_uniform"
	#iterative_prior_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_iterative_variant_gene_prior_v3_pip_level_bootstrapped.txt")
	#iterative_sampler_violin_plot2 <- make_violin_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)


	# Join plots together
	#pp <- plot_grid(iterative_pmces_violin_plot, iterative_sampler_violin_plot, iterative_sampler_violin_plot2, ncol=1)
	#output_file <- paste0(visualize_tgfm_dir, "tissue_violinplot_of_distribution_prior_probabilities_joint_", trait_name_readable,".pdf")
	#ggsave(pp, file=output_file, width=10.2, height=8.6, units="in")

	# Iterative version
	output_file <- paste0(visualize_tgfm_dir, "tissue_violinplot_of_distribution_prior_probabilities_", trait_name_readable,".pdf")
	ggsave(iterative_sampler_violin_plot, file=output_file, width=10.2, height=4.6, units="in")


	# TGLR version
	#bs_nn_sldsc_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt")
	#bs_nn_tglr_violin_plot <- make_violin_plot_showing_distribution_of_bootstrapped_taus_of_each_tissue(bs_nn_sldsc_file, "bootstrapped nn TGLR", trait_name_readable, threshold=1e-8)
	#output_file <- paste0(visualize_tgfm_dir, "tissue_violinplot_of_distribution_bs_nn_tglr_", trait_name_readable,".pdf")
	#ggsave(bs_nn_tglr_violin_plot, file=output_file, width=10.2, height=4.6, units="in")
}


if (FALSE) {


##########################################################
# Swarm-plot showing gene-tissue PIPs colored by tissue group for each trait
##########################################################
method_version="susie_sampler_uniform"
beeswarm_plot <- make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait(trait_names, trait_names_readable, method_version, tgfm_results_dir)
output_file <- paste0(visualize_tgfm_dir, "beeswarm_gene_tissue_pip_colored_by_tissue_group_", method_version,".pdf")
ggsave(beeswarm_plot, file=output_file, width=7.2, height=4.5, units="in")
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
beeswarm_plot <- make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait(trait_names, trait_names_readable, method_version, tgfm_results_dir)
output_file <- paste0(visualize_tgfm_dir, "beeswarm_gene_tissue_pip_colored_by_tissue_group_", method_version,".pdf")
ggsave(beeswarm_plot, file=output_file, width=7.2, height=4.5, units="in")
method_version="susie_sampler_tglr_bootstrapped_nonnegative_sampler"
beeswarm_plot <- make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait(trait_names, trait_names_readable, method_version, tgfm_results_dir)
output_file <- paste0(visualize_tgfm_dir, "beeswarm_gene_tissue_pip_colored_by_tissue_group_", method_version,".pdf")
ggsave(beeswarm_plot, file=output_file, width=7.2, height=4.5, units="in")


##########################################################
# Bar plot showing expected number of causal tissue categories
##########################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_tissue_categories_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_tissue_categories_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_bar_plot, file=output_file, width=7.2, height=4.5, units="in")

}


##########################################################
# Bar plot showing expected number of causal gene-tissue pairs in each tissue
##########################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_uniform"
	tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_gene_tissue_pairs_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_bar_plot, file=output_file, width=7.2, height=3.7, units="in")

	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_gene_tissue_pairs_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_bar_plot, file=output_file, width=7.2, height=3.7, units="in")

	# PMCES approach
	method_version="susie_sampler_tglr_bootstrapped_nonnegative_sampler"
	tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_gene_tissue_pairs_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_bar_plot, file=output_file, width=7.2, height=3.7, units="in")	
}


##########################################################
# Bar plot TGFM causal tissue p-values
##########################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	tissue_pvalue_bar_plot <- make_bar_plot_showing_tgfm_tissue_pvalues_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_tgfm_log_pvalue_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_pvalue_bar_plot, file=output_file, width=7.2, height=3.7, units="in")

	# Sampler approach
	method_version="susie_sampler_uniform"
	tissue_pvalue_bar_plot <- make_bar_plot_showing_tgfm_tissue_pvalues_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_tgfm_log_pvalue_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_pvalue_bar_plot, file=output_file, width=7.2, height=3.7, units="in")

}



##########################################################
# Barplot with standard errors showing fraction of mediated components across traits
##########################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_average_expression_mediated_probability_se_barplot.pdf")
med_prob_se_barplot <- make_mediated_prob_se_barplot(trait_names, trait_names_readable, method_version, tgfm_results_dir)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")

##########################################################
# Barplot with standard errors showing fraction of high pip genetic elements from gene expression
##########################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_average_fraction_of_genetic_elements_from_gene_expression_se_barplot.pdf")
med_prob_se_barplot <- make_fraction_of_genetic_elements_from_gene_expression_se_barplot(trait_names, trait_names_readable, method_version, tgfm_results_dir)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")
}


if (FALSE) {
##########################################################
# Heatmap showing Jaccard Index overlaps
##########################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"


	tissue_overlap_heatmap <- make_tissue_overlap_jaccard_index_heatmap(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tgfm_tissue_overlap_jaccard_index_heatmap_", trait_name_readable, "_", method_version, ".pdf")
	ggsave(tissue_overlap_heatmap, file=output_file, width=7.2, height=5.5, units="in")
}
}
























































if (FALSE) {
##################################################
# Bar plot showing iterative prior probabilities given distribution estimates
##################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# PMCES approach
	method_version="susie_pmces_variant_gene"
	iterative_prior_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_iterative_variant_gene_prior_bootstrapped.txt")
	barplot1 <- make_bar_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)

	# S-ldsc approach
	sldsc_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_organized_res.txt")
	barplot2 <- make_bar_plot_showing_sldsc_tau_of_each_tissue(sldsc_file, trait_name_readable)

	pp <- plot_grid(barplot1, barplot2, ncol=1)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_distribution_iterative_prior_probabilities_", trait_name_readable, "_", method_version,".pdf")
	ggsave(pp, file=output_file, width=7.2, height=8.6, units="in")

	# Sampler approach
	method_version="susie_sampler_variant_gene"
	iterative_prior_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_iterative_variant_gene_prior_bootstrapped.txt")
	barplot1 <- make_bar_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)

	# S-ldsc approach
	sldsc_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_organized_res.txt")
	barplot2 <- make_bar_plot_showing_sldsc_tau_of_each_tissue(sldsc_file, trait_name_readable)

	pp <- plot_grid(barplot1, barplot2, ncol=1)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_distribution_iterative_prior_probabilities_", trait_name_readable, "_", method_version,".pdf")
	ggsave(pp, file=output_file, width=7.2, height=8.6, units="in")


}


##################################################
# Bar plot showing iterative prior probabilities given point estimates of priors
##################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_variant_gene"
	iterative_prior_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_iterative_variant_gene_prior.txt")
	barplot1 <- make_bar_plot_showing_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)
	#output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_iterative_prior_probabilities_", trait_name_readable, "_", method_version,".pdf")
	#ggsave(barplot1, file=output_file, width=7.2, height=4.2, units="in")


	# S-ldsc approach
	sldsc_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_organized_res.txt")
	barplot2 <- make_bar_plot_showing_sldsc_tau_of_each_tissue(sldsc_file, trait_name_readable)

	pp <- plot_grid(barplot1, barplot2, ncol=1)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_iterative_prior_probabilities_", trait_name_readable, "_", method_version,".pdf")
	ggsave(pp, file=output_file, width=7.2, height=8.6, units="in")
}


}




if (FALSE) {
##########################################################
# Scatter plot comparing number of pip genes in each tissue with tissue's expression-mediated h2
##########################################################

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_variant_gene"
	scatterplot <- make_scatterplot_comparing_tglr_expression_h2_and_n_tgfm_pip_across_tissues(trait_name, trait_name_readable, method_version, tgfm_results_dir)
	output_file <- paste0(visualize_tgfm_dir, "scatter_comparing_tglr_expression_h2_and_n_tgfm_pip_across_tissues_", trait_name_readable, "_", method_version,".pdf")
	ggsave(scatterplot, file=output_file, width=7.2, height=4.5, units="in")
}

}





if (FALSE) {
# Extract trait names
trait_df <- read.table(independent_trait_names_file, header=TRUE, sep="\t")
trait_names <- as.character(trait_df$study_name)
trait_names <- c("biochemistry_Cholesterol", "biochemistry_VitaminD", "blood_MEAN_PLATELET_VOL", "blood_MONOCYTE_COUNT", "body_BMIz", "body_WHRadjBMIz", "bp_DIASTOLICadjMEDz")
trait_names_readable <- c("Cholesterol", "VitaminD", "Platelet_vol", "Monocyte_count", "BMI", "WHRadjBMI", "Diastolic_BP")



##########################################################
# Swarm-plot showing gene-tissue PIPs colored by tissue group for each trait
##########################################################
method_version="susie_pmces_sparse_variant_gene_tissue"
#beeswarm_plot <- make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait(trait_names, trait_names_readable, method_version, tgfm_results_dir)
output_file <- paste0(visualize_tgfm_dir, "beeswarm_gene_tissue_pip_colored_by_tissue_group_", method_version,".pdf")
#ggsave(beeswarm_plot, file=output_file, width=7.2, height=4.5, units="in")
method_version="susie_sampler_sparse_variant_gene_tissue"
beeswarm_plot <- make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait(trait_names, trait_names_readable, method_version, tgfm_results_dir)
output_file <- paste0(visualize_tgfm_dir, "beeswarm_gene_tissue_pip_colored_by_tissue_group_", method_version,".pdf")
ggsave(beeswarm_plot, file=output_file, width=7.2, height=4.5, units="in")



##########################################################
# Bar plot showing expected number of causal gene-tissue pairs in each tissue
##########################################################

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_sparse_variant_gene_tissue"
	tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_gene_tissue_pairs_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_bar_plot, file=output_file, width=7.2, height=4.5, units="in")

	# PMCES approach
	method_version="susie_pmces_sparse_variant_gene_tissue"
	#tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_gene_tissue_pairs_", trait_name_readable, "_", method_version,".pdf")
	#ggsave(tissue_bar_plot, file=output_file, width=7.2, height=4.5, units="in")
}


##########################################################
# Scatter plot comparing number of pip genes in each tissue with tissue's expression-mediated h2
##########################################################

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_sparse_variant_gene_tissue"
	scatterplot <- make_scatterplot_comparing_tglr_expression_h2_and_n_tgfm_pip_across_tissues(trait_name, trait_name_readable, method_version, tgfm_results_dir)
	output_file <- paste0(visualize_tgfm_dir, "scatter_comparing_tglr_expression_h2_and_n_tgfm_pip_across_tissues_", trait_name_readable, "_", method_version,".pdf")
	ggsave(scatterplot, file=output_file, width=7.2, height=4.5, units="in")
}

}



































if (FALSE) {
# Get expression mediated heritabilities
anno_size = extract_anno_size(processed_tgfm_sldsc_data_dir)
# For sparse model
methods <- c("sparse_ard_no_geno_regularization")
expr_med_frac_sparse_model_df = extract_fraction_of_h2_mediated_by_gene_expression_for_sparse_models(trait_names, tgfm_sldsc_results_dir, methods, anno_size, nonneg=FALSE)
expr_med_frac_sparse_nonneg_model_df = extract_fraction_of_h2_mediated_by_gene_expression_for_sparse_models(trait_names, tgfm_sldsc_results_dir, methods, anno_size, nonneg=TRUE)

# for full model
methods <- c("baselineLD_no_qtl")
expr_med_frac_df = extract_fraction_of_h2_mediated_by_gene_expression_for_several_methods(trait_names, tgfm_sldsc_results_dir, methods)

# Extract df containing fraction of components mediated by gene expression for each trait
fraction_component_mediated_df <- extract_data_frame_containing_fraction_of_mediated_components_per_trait(trait_names, tgfm_results_dir)

# Extract fraction of components coming from each tissue
tissue_fraction_component_mediated_df <- extract_data_frame_containing_fraction_of_mediated_components_in_each_tissue_per_trait(trait_names, tgfm_results_dir)
methods <- c("sparse_ard_no_geno_regularization")
tissue_fraction_h2_mediated_sparse_nonneg_df <- extract_data_frame_containing_fraction_of_mediated_h2_in_each_tissue_per_trait(trait_names, tgfm_sldsc_results_dir, methods, anno_size, nonneg=TRUE)


##########################################################
# Stacked barplot showing tgfm tissue-fraction mediated for each trati
##########################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	
	# TGFM Components
	trait_tissue_fraction_component_mediated_df = tissue_fraction_component_mediated_df[as.character(tissue_fraction_component_mediated_df$trait)==trait_name,]
	tgfm_stacked_barplot <- stacked_barplot_for_single_trait_breaking_down_expression_effects_per_tissue(trait_tissue_fraction_component_mediated_df, gtex_colors_df, trait_name, "Fraction of expression-mediated TGFM components")
	output_file <- paste0(visualize_tgfm_dir, "stacked_barplot_showing_tgfm_per_tissue_contribution_", trait_name, ".pdf")
	ggsave(tgfm_stacked_barplot, file=output_file, width=7.2, height=4.5, units="in")

	# h2
	trait_tissue_fraction_h2_mediated_sparse_nonneg_df = tissue_fraction_h2_mediated_sparse_nonneg_df[as.character(tissue_fraction_h2_mediated_sparse_nonneg_df$trait)==trait_name,]
	h2_stacked_barplot <- stacked_barplot_for_single_trait_breaking_down_expression_effects_per_tissue(trait_tissue_fraction_h2_mediated_sparse_nonneg_df, gtex_colors_df, trait_name, "Fraction of expression-mediated heritability")
	output_file <- paste0(visualize_tgfm_dir, "stacked_barplot_showing_ldsc_per_tissue_contribution_", trait_name, ".pdf")
	ggsave(h2_stacked_barplot, file=output_file, width=7.2, height=4.5, units="in")

	# Merge two plots together with cowplot
	joint_stacked_barplot <- plot_grid(h2_stacked_barplot, tgfm_stacked_barplot, ncol=1)
	output_file <- paste0(visualize_tgfm_dir, "stacked_barplot_showing_joint_per_tissue_contribution_", trait_name, ".pdf")
	ggsave(joint_stacked_barplot, file=output_file, width=7.2, height=5.5, units="in")



}



##########################################################
# Histogram showing distribution of mediated probabilities across components (seperately for each trait)
##########################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	expression_med_prob <- extract_expression_mediated_probabilities_across_components(trait_name, tgfm_results_dir)

	histo <- plot_distribution_of_tgfm_expression_mediated_probabilities(expression_med_prob, trait_name)
	output_file <- paste0(visualize_tgfm_dir, "distribution_of_mediated_probabilities_", trait_name, ".pdf")
	ggsave(histo, file=output_file, width=7.2, height=4.5, units="in")
}


##########################################################
# Scatterplot showing correlation in expression-mediated heritability between full and sparse model
##########################################################
scatter <- make_scatterplot_comparing_expression_mediated_h2_across_traits(expr_med_frac_df, expr_med_frac_sparse_model_df, "expression-mediated h2\nfull-model", "expression-mediated h2\nsparse model")
output_file <- paste0(visualize_tgfm_dir, "expression_mediated_h2_scatter_across_traits_for_full_and_sparse_model.pdf")
ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")

scatter <- make_scatterplot_comparing_expression_mediated_h2_across_traits(expr_med_frac_sparse_nonneg_model_df, expr_med_frac_sparse_model_df, "Expression-mediated h2\nsparse nonneg model", "expression-mediated h2\nsparse model")
output_file <- paste0(visualize_tgfm_dir, "expression_mediated_h2_scatter_across_traits_for_sparse_and_sparse_nonneg_model.pdf")
ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")

##########################################################
# Scatterplot showing correlation in expression-mediated heritability and average mediated probability
##########################################################
scatter <- make_scatterplot_comparing_expression_mediated_h2_and_average_component_mediated_probability_across_traits(expr_med_frac_sparse_model_df, fraction_component_mediated_df, "Expression-mediated h2\nsparse model", "Average expression-mediated probability\nacross TGFM components")
output_file <- paste0(visualize_tgfm_dir, "expression_mediated_h2_vs_tgfm_average_mediated_prob_scatter_across_traits.pdf")
ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")

scatter <- make_scatterplot_comparing_expression_mediated_h2_and_average_component_mediated_probability_across_traits(expr_med_frac_sparse_nonneg_model_df, fraction_component_mediated_df, "Expression-mediated h2\nsparse nonneg model", "Average expression-mediated probability\nacross TGFM components")
output_file <- paste0(visualize_tgfm_dir, "expression_mediated_nonneg_h2_vs_tgfm_average_mediated_prob_scatter_across_traits.pdf")
ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")


}


