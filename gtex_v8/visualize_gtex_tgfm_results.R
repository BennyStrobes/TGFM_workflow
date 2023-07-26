args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(reshape2)
library(ggbeeswarm)
options(warn=1)
options(bitmapType='cairo')


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

make_number_of_high_pip_elements_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, pip_threshold, independent_traits) {
	trait_names_vec <- c()
	number_elements_vec <- c()
	element_type_vec <- c()
	mediated_fraction <- c()

	variant_name = paste0("Variant (PIP >= ", pip_threshold, ")")
	gene_name = paste0("Gene-Tissue (PIP >= ", pip_threshold, ")")


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {


		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_threshold, ".txt")
		tmp <- read.table(summary_file, header=TRUE)


		variant_count = tmp$count[1]
		gene_count = tmp$count[2]

		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		number_elements_vec <- c(number_elements_vec, variant_count)
		number_elements_vec <- c(number_elements_vec, gene_count)
		element_type_vec <- c(element_type_vec, variant_name)
		element_type_vec <- c(element_type_vec, gene_name)
		mediated_fraction <- c(mediated_fraction, gene_count/(gene_count+variant_count))
		}
	}

	df <- data.frame(trait=trait_names_vec, number_of_elements=number_elements_vec, element_type=factor(element_type_vec, levels=c(gene_name, variant_name)))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$element_type) == variant_name,]

	indices = order(-tmp_df$number_of_elements)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=number_of_elements, fill=element_type)) +
    		#geom_bar(stat="identity", position=position_dodge()) +\
    		geom_bar(stat="identity", position="stack") +
    		scale_fill_manual(values=c("dodgerblue3", "grey")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="No. fine-mapped\ngenetic elements", x="", fill="") +
    		figure_theme() +
    		theme(legend.position=c(0.8, 0.86))
    return(p)
}

make_number_of_high_pip_cross_2_threshold_gene_tissue_pairs_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, pip_thresholds, independent_traits) {
	trait_names_vec <- c()
	number_elements_vec <- c()
	element_type_vec <- c()
	mediated_fraction <- c()

	threshold1_name <- paste0(pip_thresholds[1], " <= PIP < ", pip_thresholds[2])
	threshold2_name <- paste0("PIP  >= ", pip_thresholds[2])


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {


		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_thresholds[1], ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		variant_count1 = tmp$count[1]
		gene_count1 = tmp$count[2]

		mediated_fraction <- c(mediated_fraction, gene_count1/(gene_count1+variant_count1))

		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_thresholds[2], ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		variant_count2 = tmp$count[1]
		gene_count2 = tmp$count[2]


		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		number_elements_vec <- c(number_elements_vec, gene_count1 - gene_count2)
		number_elements_vec <- c(number_elements_vec, gene_count2)
		element_type_vec <- c(element_type_vec, threshold1_name)
		element_type_vec <- c(element_type_vec, threshold2_name)
		}
	}

	df <- data.frame(trait=trait_names_vec, number_of_elements=number_elements_vec, element_type=factor(element_type_vec, levels=c(threshold1_name, threshold2_name)))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$element_type) == threshold1_name,]

	indices = order(-tmp_df$number_of_elements)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=number_of_elements, fill=element_type)) +
    		#geom_bar(stat="identity", position=position_dodge()) +\
    		geom_bar(stat="identity", position="stack") +
    		#scale_fill_manual(values=c("dodgerblue3", "grey")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="No. fine-mapped\ngenetic elements", x="", fill="") +
    		figure_theme() +
    		theme(legend.position=c(0.8, 0.86))
    return(p)
}

make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2 <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits) {
	trait_names_vec <- c()
	number_elements_vec <- c()
	pip_threshold_vec <- c()

	total_hits_vec <- c()
	mid_point_vec <- c()
	trait_names_vec2 <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_cross_pip_threshold.txt")
			tmp_data <- read.table(summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="gene",]

			tmp_data = tmp_data[tmp_data$PIP_threshold >= .2,]
			
			number_elements_vec <- c(number_elements_vec, tmp_data$n_elements)
			pip_threshold_vec <- c(pip_threshold_vec, tmp_data$PIP_threshold)
			trait_names_vec <- c(trait_names_vec, rep(trait_name_readable, length(tmp_data$PIP_threshold)))

			total_hits_vec <- c(total_hits_vec, sum(tmp_data$n_elements))
			trait_names_vec2 <- c(trait_names_vec2, trait_name_readable)
			mid_point_vec <- c(mid_point_vec, sum(tmp_data$n_elements[tmp_data$PIP_threshold >= .5]))

		}

	}
	df <- data.frame(trait=trait_names_vec, PIP=pip_threshold_vec, num_elements=number_elements_vec)

	indices = order(-total_hits_vec)
	df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])


	df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)

   p <- ggplot() +
  	geom_bar(data=df, aes(x = trait,y = num_elements,fill = PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	scale_y_continuous(breaks=c(0.0,sqrt(5), sqrt(20), sqrt(50), sqrt(100), sqrt(200), sqrt(400), sqrt(600)), labels=c(0, 5, 20, 50, 100, 200, 400,600)) +
  	#theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
  	theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	#scale_fill_gradient(low="lightsteelblue3", high="dodgerblue3") +
  	scale_fill_viridis_c(option = "C") +
  	 geom_errorbar(data  = df2, aes(x=trait,y=midpoint, ymax=midpoint, ymin=midpoint)) +
  	 labs(x="", y="No. fine-mapped\ngene-tissue pairs")
 	return(p)
}


make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits) {
	trait_names_vec <- c()
	number_elements_vec <- c()
	pip_threshold_vec <- c()

	total_hits_vec <- c()
	mid_point_vec <- c()
	trait_names_vec2 <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_cross_pip_threshold.txt")
			tmp_data <- read.table(summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="gene",]
			
			number_elements_vec <- c(number_elements_vec, tmp_data$n_elements)
			pip_threshold_vec <- c(pip_threshold_vec, tmp_data$PIP_threshold)
			trait_names_vec <- c(trait_names_vec, rep(trait_name_readable, length(tmp_data$PIP_threshold)))

			total_hits_vec <- c(total_hits_vec, sum(tmp_data$n_elements))
			trait_names_vec2 <- c(trait_names_vec2, trait_name_readable)
			mid_point_vec <- c(mid_point_vec, sum(tmp_data$n_elements[tmp_data$PIP_threshold >= .5]))

		}

	}
	df <- data.frame(trait=trait_names_vec, PIP=pip_threshold_vec, num_elements=number_elements_vec)

	indices = order(-total_hits_vec)
	df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])


	df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)

   p <- ggplot() +
  	geom_bar(data=df, aes(x = trait,y = num_elements,fill = PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	scale_y_continuous(breaks=c(0.0,sqrt(5), sqrt(20), sqrt(50), sqrt(100), sqrt(200), sqrt(400), sqrt(600)), labels=c(0, 5, 20, 50, 100, 200, 400,600)) +
  	#theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
  	theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	#scale_fill_gradient(low="lightsteelblue3", high="dodgerblue3") +
  	scale_fill_viridis_c(option = "C") +
  	 geom_errorbar(data  = df2, aes(x=trait,y=midpoint, ymax=midpoint, ymin=midpoint)) +
  	 labs(x="", y="No. fine-mapped\ngene-tissue pairs")
 	return(p)
}


make_number_of_high_pip_cross_3_threshold_gene_tissue_pairs_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, pip_thresholds, independent_traits) {
	trait_names_vec <- c()
	number_elements_vec <- c()
	element_type_vec <- c()
	
	mediated_fraction <- c()

	trait_names_vec2 <- c()
	number_elements_vec2 <- c()

	#threshold1_name <- paste0(pip_thresholds[1], " <= PIP < ", pip_thresholds[2])
	#threshold2_name <- paste0(pip_thresholds[2], " <= PIP < ", pip_thresholds[3])
	#threshold3_name <- paste0("PIP  >= ", pip_thresholds[3])
	threshold1_name <- paste0("[ ", pip_thresholds[1],  ", ", pip_thresholds[2], ")")
	threshold2_name <- paste0("[ ", pip_thresholds[2],  ", ", pip_thresholds[3], ")")
	threshold3_name <- paste0("[ ", pip_thresholds[3],  ", 1.0]")


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {


		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_thresholds[1], ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		variant_count1 = tmp$count[1]
		gene_count1 = tmp$count[2]

		mediated_fraction <- c(mediated_fraction, gene_count1/(gene_count1+variant_count1))

		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_thresholds[2], ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		variant_count2 = tmp$count[1]
		gene_count2 = tmp$count[2]

		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_thresholds[3], ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		variant_count3 = tmp$count[1]
		gene_count3 = tmp$count[2]

		trait_names_vec2 <- c(trait_names_vec2, trait_name_readable)
		number_elements_vec2 <- c(number_elements_vec2, gene_count1)

		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		number_elements_vec <- c(number_elements_vec, gene_count1 - gene_count2)
		number_elements_vec <- c(number_elements_vec, gene_count2 - gene_count3)
		number_elements_vec <- c(number_elements_vec, gene_count3)
		element_type_vec <- c(element_type_vec, threshold1_name)
		element_type_vec <- c(element_type_vec, threshold2_name)
		element_type_vec <- c(element_type_vec, threshold3_name)
		}
	}

	df <- data.frame(trait=trait_names_vec, number_of_elements=number_elements_vec, element_type=factor(element_type_vec, levels=c(threshold3_name,threshold2_name, threshold1_name)))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = data.frame(trait=trait_names_vec2, number_of_elements=number_elements_vec2)

	indices = order(-tmp_df$number_of_elements)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	new_mediated_fraction = round(mediated_fraction*100, digits = 0)
	mediated_fraction_string <- c()
	for (ii in 1:length(new_mediated_fraction)) {
		mediated_fraction_string <- c(mediated_fraction_string, paste0(new_mediated_fraction[ii], "%"))
	}


	df2 <- data.frame(traity=trait_names_vec2, number_of_elementsy=number_elements_vec2, fractiony=mediated_fraction_string)
	df2$trait <- recode(df2$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")
	df2$trait = factor(df2$trait, levels=as.character(tmp_df$trait)[indices])


	maxy = max(number_elements_vec2) + 2.5


	p <- ggplot() +
    		#geom_bar(stat="identity", position=position_dodge()) +\
    		geom_bar(data=df, aes(x=trait, y=number_of_elements, fill=element_type), stat="identity", position="stack") +
    		#scale_fill_manual(values=c("dodgerblue3", "grey")) +
    		scale_fill_manual(values=c("dodgerblue3", "steelblue1", "lightsteelblue3"))+
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		#geom_text(data=df2, aes(x = traity, y = number_of_elementsy, label = fractiony), position = position_stack(vjust = 0.5))+
    		geom_text(data=df2, aes(x = traity, y = number_of_elementsy, label = fractiony), vjust=-0.1)+
    		labs(y="No. fine-mapped\ngenetic elements", x="", fill="Gene-Tissue PIP") +
    		figure_theme() +
    		theme(legend.position=c(0.89, 0.73)) +
    		scale_y_continuous(limits = c(0, maxy))
    return(p)
}


make_number_of_high_pip_gene_tissue_pairs_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, pip_threshold, independent_traits) {
	trait_names_vec <- c()
	number_elements_vec <- c()
	element_type_vec <- c()
	mediated_fraction <- c()

	gene_name = paste0("Gene-Tissue (PIP >= ", pip_threshold, ")")


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {


		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_threshold, ".txt")
		tmp <- read.table(summary_file, header=TRUE)


		variant_count = tmp$count[1]
		gene_count = tmp$count[2]

		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		number_elements_vec <- c(number_elements_vec, gene_count)
		element_type_vec <- c(element_type_vec, gene_name)
		mediated_fraction <- c(mediated_fraction, gene_count/(gene_count+variant_count))
		}
	}

	df <- data.frame(trait=trait_names_vec, number_of_elements=number_elements_vec, element_type=factor(element_type_vec))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$element_type) == gene_name,]

	indices = order(-tmp_df$number_of_elements)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=number_of_elements)) +
    		#geom_bar(stat="identity", position=position_dodge()) +\
    		geom_bar(stat="identity", position="stack") +
    		#scale_fill_manual(values=c("dodgerblue3", "grey")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="No. fine-mapped\ngenetic elements", x="", fill="") +
    		figure_theme() +
    		theme(legend.position=c(0.8, 0.86))
    return(p)
}

make_expected_fraction_of_genetic_elements_from_gene_expression_se_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, independent_traits) {

	trait_names_vec <- c()
	fraction_mediated_vec <- c()
	fraction_mediated_se_vec <- c()



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {

			summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_expected_n_causal_genetic_elements.txt")
			tmp <- read.table(summary_file, header=TRUE)


			variant_count = tmp$count[1]
			gene_count = tmp$count[2]
			total_count = gene_count + variant_count

			mean_mediated_prob <- gene_count/(total_count)
		
			se_mediated_prob = sqrt((mean_mediated_prob*(1.0-mean_mediated_prob))/(total_count))
			trait_names_vec <- c(trait_names_vec, trait_name_readable)
			fraction_mediated_vec <- c(fraction_mediated_vec, mean_mediated_prob)
			fraction_mediated_se_vec <- c(fraction_mediated_se_vec, se_mediated_prob)
		}
	}

	df <- data.frame(trait=trait_names_vec, average_mediated_probability=fraction_mediated_vec, average_mediated_probability_se=fraction_mediated_se_vec)
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	indices = order(df$average_mediated_probability)
	df$trait = factor(df$trait, levels=as.character(df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=average_mediated_probability)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		#scale_fill_manual(values=c("plum2", "orchid2", "magenta2")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="Expected fraction of PIP\nmediated by gene expression", x="") +
    		geom_errorbar( aes(ymin=average_mediated_probability-(1.96*average_mediated_probability_se), ymax=average_mediated_probability+(1.96*average_mediated_probability_se)), width=0.2, colour="grey50", position=position_dodge(.9)) +
    		figure_theme() +
    		theme(legend.position="bottom")
    return(p)
}


make_fraction_of_genetic_elements_from_gene_expression_se_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, independent_traits) {

	trait_names_vec <- c()
	fraction_mediated_vec <- c()
	fraction_mediated_se_vec <- c()
	pip_threshold_vec <- c()

	pip_thresholds <- c("0.25", "0.5", "0.75")


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {

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
	}

	df <- data.frame(trait=trait_names_vec, average_mediated_probability=fraction_mediated_vec, average_mediated_probability_se=fraction_mediated_se_vec, pip=factor(pip_threshold_vec,levels=pip_thresholds))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$pip) == "0.25",]

	indices = order(tmp_df$average_mediated_probability)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=average_mediated_probability, fill=pip)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		scale_fill_manual(values=c("plum2", "orchid2", "magenta2")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="Fraction of fine-mapped\ngenetic elements mediated\nby gene expression", x="") +
    		geom_errorbar( aes(ymin=average_mediated_probability-(1.96*average_mediated_probability_se), ymax=average_mediated_probability+(1.96*average_mediated_probability_se)), width=0.2, colour="grey50", position=position_dodge(.9)) +
    		figure_theme() +
    		theme(legend.position="bottom")
    return(p)
}
make_tissue_mediated_prob_se_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, lb_thresh=.1) {

	trait_names_vec <- c()
	fraction_mediated_vec <- c()
	fraction_mediated_se_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]

		print(trait_name_readable)
		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_component_level_summary.txt")
		
		tmp <- read.table(summary_file, header=TRUE)
		mediated_probs <- 1.0 - tmp$non_mediated_probability

		tmp = tmp[,5:dim(tmp)[2]]
		print(sort(colMeans(tmp))/sum(colMeans(tmp)))
		
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
    		geom_bar( aes(x=trait, y=average_mediated_probability), stat="identity", fill="skyblue", alpha=0.7, width=.75) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    		labs(y="Average expression\nmediated probability", x="") +
    		geom_errorbar( aes(x=trait, ymin=average_mediated_probability-(1.96*average_mediated_probability_se), ymax=average_mediated_probability+(1.96*average_mediated_probability_se)), width=0.4, colour="orange", alpha=0.9, size=.9) +
    		figure_theme()
    return(p)
}

make_mediated_prob_se_scatterplot <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, tgfm_sldsc_results_dir) {

	trait_names_vec <- c()
	fraction_mediated_vec <- c()
	fraction_mediated_se_vec <- c()
	h2_mediated_vec <- c()
	h2_mediated_se_vec <- c()
	if (FALSE) {
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name != "CAD") {

		summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_component_level_summary.txt")
		
		tmp <- read.table(summary_file, header=TRUE)
		mediated_probs <- 1.0 - tmp$non_mediated_probability
		
		mean_mediated_prob = mean(mediated_probs)
		se_mediated_prob = sqrt((mean_mediated_prob*(1.0-mean_mediated_prob))/(length(mediated_probs)/100.0))
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		fraction_mediated_vec <- c(fraction_mediated_vec, mean_mediated_prob)
		fraction_mediated_se_vec <- c(fraction_mediated_se_vec, se_mediated_prob)


		h2_summary_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_h2_5_50_med.txt")
		tmp <- read.table(h2_summary_file, header=TRUE, sep="\t")
		h2_mediated_vec <- c(h2_mediated_vec, tmp$h2_med[1])
		h2_mediated_se_vec <- c(h2_mediated_se_vec, tmp$h2_med_se[1])
		}

	}
	df <- data.frame(trait=trait_names_vec, average_mediated_probability=fraction_mediated_vec, average_mediated_probability_se=fraction_mediated_se_vec, mediated_h2=h2_mediated_vec, mediated_h2_se=h2_mediated_se_vec)
	}
	load(file="df.Rda")

	p <- ggplot(df, aes(x=average_mediated_probability, y=mediated_h2)) + geom_point() + figure_theme() +geom_abline(slope=1)

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
	#save(df,file="data.Rda")
	#print("DONE")
	#load("data.Rda")


	#df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	indices = order(df$average_mediated_probability)
	df$trait = factor(df$trait, levels=as.character(df$trait)[indices])

	p <- ggplot(df) +
    		geom_bar( aes(x=trait, y=average_mediated_probability), stat="identity", fill="skyblue", alpha=0.7, width=.75) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    		labs(y="Expression\nmediated prob", x="") +
    		geom_errorbar( aes(x=trait, ymin=average_mediated_probability-(1.96*average_mediated_probability_se), ymax=average_mediated_probability+(1.96*average_mediated_probability_se)), width=0.4, colour="orange", alpha=0.9, size=.9) +
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

	df <- df[2:nrow,]
	df$tissue = df$element_name

	tissue_vec <- c()
	prior_prob_vec <- c()
	siggy <- c()

	for (tissue_iter in 1:length(df$tissue)) {
		tissue_name <- as.character(df$tissue[tissue_iter])
		prob_string = as.character(df$prior_distribution[tissue_iter])
		prob_vec = as.numeric(strsplit(prob_string,";")[[1]])

		zz = mean(prob_vec)/sd(prob_vec)

		#pvalue = 2*pnorm(q=zz, lower.tail=FALSE)
		#print(pvalue)
		pvalue = sum(prob_vec < threshold)/length(prob_vec)
		#print(pvalue2)
		#print(paste0(pvalue, "   ", pvalue2))
		if (is.na(pvalue)) {
			pvalue = 1.0
		}

		if (pvalue < (.1/37)) {
			sigger = 'significant'
		} else {
			sigger = 'non_significant'
		}


		prior_prob_vec <- c(prior_prob_vec, prob_vec)
		tissue_vec <- c(tissue_vec, rep(tissue_name, length(prob_vec)))
		siggy <- c(siggy, rep(sigger, length(prob_vec)))
	}


	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	tissue_vec = str_replace_all(as.character(tissue_vec), "-", "_")


	df2 <- data.frame(tissue=factor(tissue_vec, levels=as.character(df$tissue)), probability=prior_prob_vec, significance=factor(siggy))
	#df2$tissue = str_replace_all(as.character(df2$tissue), "-", "_")
	df2$tissue <- recode(df2$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	#df2$tissue = factor(df2$tissue)

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


make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh, selected_traits) {
	tissue_vec <- c()
	count_vec <- c()
	trait_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		if (trait_name %in% selected_traits) {
		trait_name_readable <- trait_names_readable[trait_iter]
		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")


		tmper <- c()
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		if (sum(trait_df$PIP > .1) > 0) {
			#print(head(tmp_df))
			for (tissue_iter in 1:length(ordered_tissues)) {
				tissue_name <- ordered_tissues[tissue_iter]
				pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				#count_vec <- c(count_vec, pip_sum)
				trait_vec <- c(trait_vec, trait_name_readable)
				tmper <- c(tmper, pip_sum)
			}
			tmper = tmper/sum(tmper)
			for (itera in 1:length(tmper)) {
				tissue_name <- ordered_tissues[itera]
				count_vec <- c(count_vec, tmper[itera])
			}
		}
		}
	}

	new_trait_names = sort(unique(trait_vec))

	valid_tissues <- c()
	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues), trait=factor(trait_vec, levels=new_trait_names))
	for (tissue_iter in 1:length(ordered_tissues)) {
		tissue_name = ordered_tissues[tissue_iter]
		indices = as.character(df$tissue) == tissue_name
		if (max(df$expected_causal_genes[indices]) > .18) {
			valid_tissues <- c(valid_tissues, tissue_name)
		}
	}


	df = df[as.character(df$tissue) %in% valid_tissues,]

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin (sun exposed)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord", Whole_Blood="Whole Blood", Colon_Sigmoid="Colon Sigmoid", Artery_Tibial="Artery Tibial", Brain_BasalGanglia="Brain Basal Ganglia", Esophagus_Mucosa="Esophagus Mucosa", Brain_Cortex="Brain Cortex", Artery_Heart="Artery Heart", Adrenal_Gland="Adrenal Gland")
	df$tissue = factor(df$tissue)
	df$value = df$expected_causal_genes




	matrix <- dcast(df, tissue ~ trait)
	tissue_names <- matrix$tissue
	matrix <- as.matrix(matrix[,2:dim(matrix)[2]])


	#ord <- hclust( dist(matrix, method = "euclidean"), method = "ward.D" )$order
	#df$tissue <- factor(df$tissue, levels=as.character(tissue_names)[ord])
	custom_tissue_names <- c("Spleen", "Whole Blood", "Lymphocytes", "Thyroid", "Colon Sigmoid", "Pancreas", "Fibroblast", "Artery Tibial", "Lung", "Brain Basal Ganglia", "Brain Cortex", "Esophagus Mucosa", "Skin (sun exposed)", "Liver", "Artery Heart", "Adrenal Gland")
	df$tissue <- factor(df$tissue, levels=custom_tissue_names)
	
	#ord2 <- hclust( dist(t(matrix), method = "euclidean"), method = "ward.D" )$order
	#df$trait <- factor(df$trait, levels=as.character(new_trait_names)[ord2])
	custom_orderd_traits = rev(c("Autoimmune", "Monocyte count", "Hypothyroidism", "HbA1c", "Height", "BMI", "Neuroticism", "Eczema", "VitaminD", "LDL direct", "Hypertension"))
	df$trait <- factor(df$trait, levels=custom_orderd_traits)


	pp <- ggplot(df, aes(tissue, trait, fill= expected_causal_genes)) + 
  		geom_tile() +
  		figure_theme() +
  		#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low = "white", high = "dodgerblue3") +
  		labs(x="",y="", fill="Expected fraction of fine-mapped  \ngene-tissue pairs")

  	return(pp)	
}

make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits_for_figure2 <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh, valid_tissues) {
	tissue_vec <- c()
	count_vec <- c()
	trait_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")


		tmper <- c()
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		if (sum(trait_df$PIP > .7) > 0) {
			#print(head(tmp_df))
			for (tissue_iter in 1:length(ordered_tissues)) {
				tissue_name <- ordered_tissues[tissue_iter]
				pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				#count_vec <- c(count_vec, pip_sum)
				trait_vec <- c(trait_vec, trait_name_readable)
				tmper <- c(tmper, pip_sum)
			}
			tmper = tmper/sum(tmper)
			for (itera in 1:length(tmper)) {
				tissue_name <- ordered_tissues[itera]
				count_vec <- c(count_vec, tmper[itera])
			}
		}
	}

	new_trait_names = sort(unique(trait_vec))

	valid_tissues <- c()
	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues), trait=factor(trait_vec, levels=new_trait_names))
	for (tissue_iter in 1:length(ordered_tissues)) {
		tissue_name = ordered_tissues[tissue_iter]
		indices = as.character(df$tissue) == tissue_name
		if (max(df$expected_causal_genes[indices]) > .2) {
			valid_tissues <- c(valid_tissues, tissue_name)
		}
	}


	df = df[as.character(df$tissue) %in% valid_tissues,]

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$tissue = factor(df$tissue)
	df$value = df$expected_causal_genes




	matrix <- dcast(df, tissue ~ trait)
	tissue_names <- matrix$tissue
	matrix <- as.matrix(matrix[,2:dim(matrix)[2]])

	ord <- hclust( dist(matrix, method = "euclidean"), method = "ward.D" )$order
	df$tissue <- factor(df$tissue, levels=as.character(tissue_names)[ord])
	ord2 <- hclust( dist(t(matrix), method = "euclidean"), method = "ward.D" )$order

	df$trait <- factor(df$trait, levels=as.character(new_trait_names)[ord2])

	pp <- ggplot(df, aes(tissue, trait, fill= expected_causal_genes)) + 
  		geom_tile() +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  		theme(legend.position="right") +
  		scale_fill_gradient(low = "grey", high = "red") +
  		labs(x="",y="", fill="Tissue\nfraction") 

  	return(pp)


}


make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh, valid_tissues) {
	tissue_vec <- c()
	count_vec <- c()
	trait_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")


		tmper <- c()
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		if (sum(trait_df$PIP > .7) > 0) {
			#print(head(tmp_df))
			for (tissue_iter in 1:length(ordered_tissues)) {
				tissue_name <- ordered_tissues[tissue_iter]
				pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				#count_vec <- c(count_vec, pip_sum)
				trait_vec <- c(trait_vec, trait_name_readable)
				tmper <- c(tmper, pip_sum)
			}
			tmper = tmper/sum(tmper)
			for (itera in 1:length(tmper)) {
				tissue_name <- ordered_tissues[itera]
				count_vec <- c(count_vec, tmper[itera])
			}
		}
	}

	new_trait_names = sort(unique(trait_vec))

	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues), trait=factor(trait_vec, levels=new_trait_names))
	for (tissue_iter in 1:length(ordered_tissues)) {
		tissue_name = ordered_tissues[tissue_iter]
		indices = as.character(df$tissue) == tissue_name
	}



	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$tissue = factor(df$tissue)
	df$value = df$expected_causal_genes


	df = df[as.character(df$tissue) %in% valid_tissues,]



	matrix <- dcast(df, tissue ~ trait)
	tissue_names <- matrix$tissue
	matrix <- as.matrix(matrix[,2:dim(matrix)[2]])

	ord <- hclust( dist(matrix, method = "euclidean"), method = "ward.D" )$order
	df$tissue <- factor(df$tissue, levels=as.character(tissue_names)[ord])
	ord2 <- hclust( dist(t(matrix), method = "euclidean"), method = "ward.D" )$order

	df$trait <- factor(df$trait, levels=as.character(new_trait_names)[ord2])

	pp <- ggplot(df, aes(tissue, trait, fill= expected_causal_genes)) + 
  		geom_tile() +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low = "grey", high = "red")

  	return(pp)


}

make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait_for_figure <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir, ordered_tissues) {
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


	trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
	trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")


	pip_thresh = .1
	valid_tissues <- c()
	pips <- c()
	tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
	for (tissue_iter in 1:length(ordered_tissues)) {
		tissue_name <- ordered_tissues[tissue_iter]
		pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
		pips <- c(pips, pip_sum)

	}
	tmp=ordered_tissues[order(-pips)]

	valid_tissues = tmp[1:8]
	valid_tissues = unique(valid_tissues)

	pip_threshs <- c(.1, .3, .5)
	#pip_threshs <- c(.3, .5, .7)


	#ordered_tissues <- sort(unique(as.character(trait_df$tissue_name)))
	new_ordered_tissues <- c()
	tmper <- c()
	for (pip_iter in 1:length(pip_threshs)) {
		pip_thresh <- pip_threshs[pip_iter]
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		#print(head(tmp_df))
		for (tissue_iter in 1:length(ordered_tissues)) {
			tissue_name <- ordered_tissues[tissue_iter]
			if (tissue_name %in% valid_tissues) {
				new_ordered_tissues <- c(new_ordered_tissues, tissue_name)
				pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				count_vec <- c(count_vec, pip_sum)
				pip_thresh_vec <- c(pip_thresh_vec, pip_thresh)
				if (pip_thresh == pip_threshs[1]) {
					tmper <- c(tmper, pip_sum)
				}
			}
		}


	}


	ord <- order(tmper)
	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=new_ordered_tissues[ord]), pip_thresh=factor(pip_thresh_vec))


	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose Sub", Adipose_Visceral_Omentum="Adipose Visceral", Breast_Mammary_Tissue="Breast Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart Atrial",Skin_Sun_Exposed_Lower_leg="Skin-Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin-No Sun", Small_Intestine_Terminal_Ileum="Small Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain anterior cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain basal ganglia", Esophagus_Gastroesophageal_Junction="Esophagus gastro jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain Spinal cord", Whole_Blood="Whole Blood", Colon_Transverse="Colon Transverse", Brain_Cerebellum="Brain Cerebellum", Colon_Sigmoid="Colon Sigmoid", Artery_Heart="Artery Heart", Adrenal_Gland="Adrenal Gland", Muscle_Skeletal="Muscle Skeletal")
	ordered_tissues2 <- as.character(df$tissue)[1:length(unique(df$tissue))]
	df$tissue = factor(df$tissue, levels=ordered_tissues2[ord])

	p<-ggplot(df, aes(x=expected_causal_genes, y=tissue, fill=pip_thresh)) +
  		geom_bar(stat="identity",position=position_dodge())+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10)) +
  		theme(legend.position="bottom") +
  		labs(x="Expected # causal genes", y="", fill="PIP threshold",title=trait_name_readable)

  	return(p)
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
	#pip_threshs <- c(.3, .5, .7)

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
		trait_df <- trait_df[trait_df$PIP > .5, ]

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
  		geom_beeswarm(cex = .9) +
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

get_sig_tissues_tgfm_tissue_pvalues_for_single_trait <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names) {
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

	tissues = as.character(df$tissue_name[df$significant])
  	return(tissues)
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

make_enrichment_boxplot <- function(df, titler) {
	cts = unique(as.character(df$track_cell_type))
	mean_values <- c()
	for (ct_iter in 1:length(cts)) {
		ct = cts[ct_iter]
		mean_value = mean(df$orat[as.character(df$track_cell_type)==ct])
		mean_values <- c(mean_values, mean_value)
	}

	ct_ordering = order(mean_values)

	df$track_cell_type = factor(df$track_cell_type, levels=cts[ct_ordering])
	p <- ggplot(df, aes(x=track_cell_type, y=orat)) + 
  		geom_boxplot() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		labs(title=titler, x="Cell type", "Odds ratio")
  	return(p)
}

make_n_genes_vs_sample_size_scatterplot <- function(merged_tissue_cell_type_file) {
	df <- read.table(merged_tissue_cell_type_file, header=TRUE, sep="\t")
	p <- ggplot(df, aes(x=context_sample_size, y=n_gene_model_genes, color=data_set)) + geom_point() +
	figure_theme() + 
	labs(x="eQTL sample size", y="Number of genes with gene model",color="")
	return(p)
}

make_n_genes_vs_avg_n_cells_per_indi_scatterplot <- function(merged_tissue_cell_type_file) {
	df <- read.table(merged_tissue_cell_type_file, header=TRUE, sep="\t")
	df <- df[as.character(df$data_set) == "Perez_sc_pbmc",]

	p <- ggplot(df, aes(x=avg_n_cells_per_individual, y=n_gene_model_genes)) + geom_point() +
	figure_theme() + 
	labs(x="Avg. number of cells / individual", y="Number of genes with gene model",color="")
	return(p)
}



##################################
# Extract command line arguments
##################################

independent_trait_names_file <- args[1]
tgfm_sldsc_results_dir <- args[2]
tgfm_results_dir <- args[3]
tgfm_organized_results_dir <- args[4]
processed_tgfm_sldsc_data_dir <- args[5]
gtex_tissue_colors_file <- args[6]
iterative_tgfm_prior_dir <- args[7]
visualize_tgfm_dir <- args[8]

print(independent_trait_names_file)




##########################################################
# Load in data
##########################################################
# Independent traits used in https://www.nature.com/articles/s41588-020-00735-5
independent_traits <- c("body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled")



# Load in gtex tissue colors 
gtex_colors_df <- read.table(gtex_tissue_colors_file, header=TRUE, sep="\t")
gtex_colors_df$tissue_site_detail_id = as.character(gtex_colors_df$tissue_site_detail_id)
gtex_colors_df$tissue_site_detail_id[23] = "Cells_Cultured_fibroblasts"


# Extract trait names
trait_df <- read.table(independent_trait_names_file, header=TRUE, sep="\t")
trait_names <- as.character(trait_df$study_name)
trait_names_readable <- as.character(trait_df$study_name_readable)
#trait_names <- c("biochemistry_Cholesterol", "biochemistry_VitaminD", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "blood_MEAN_PLATELET_VOL", "blood_MONOCYTE_COUNT", "body_BMIz", "body_WHRadjBMIz", "bp_DIASTOLICadjMEDz", "lung_FEV1FVCzSMOKE")
#trait_names_readable <- c("Cholesterol", "VitaminD", "Reticulocyte_count", "Platelet_vol", "Monocyte_count", "BMI", "WHRadjBMI", "Diastolic_BP", "FEV1FVC")






# Extract tissue names
iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_names[1], "_component_gene_", "susie_pmces_uniform", "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
aa = read.table(iterative_prior_file, header=TRUE,sep="\t")
tissue_names = as.character(aa$element_name[2:length(aa$element_name)])

if (FALSE) {
##################################################
# Violin plot showing bootstrapped prior distributions 
##################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]
	#print(trait_name)

	# Sampler approach
	method_version="susie_pmces_uniform"
	iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
	iterative_sampler_violin_plot <- make_violin_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)

	# Iterative version
	output_file <- paste0(visualize_tgfm_dir, "tissue_violinplot_of_distribution_prior_probabilities_", trait_name_readable,".pdf")
	ggsave(iterative_sampler_violin_plot, file=output_file, width=10.2, height=4.6, units="in")


	if (trait_name != "CAD") {
	if (trait_name != "repro_ChildLess") {
	# TGLR version
	bs_nn_sldsc_file <- paste0(tgfm_sldsc_results_dir, trait_name, "_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt")
	#bs_nn_tglr_violin_plot <- make_violin_plot_showing_distribution_of_bootstrapped_taus_of_each_tissue(bs_nn_sldsc_file, "bootstrapped nn TGLR", trait_name_readable, threshold=1e-8)
	#output_file <- paste0(visualize_tgfm_dir, "tissue_violinplot_of_distribution_bs_nn_tglr_", trait_name_readable,".pdf")
	#ggsave(bs_nn_tglr_violin_plot, file=output_file, width=10.2, height=4.6, units="in")
	}
	}
}
}

##########################################################
# Bar plot TGFM causal tissue p-values
##########################################################
if (FALSE) {
valid_tissues <- c()
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	tissue_pvalue_bar_plot <- make_bar_plot_showing_tgfm_tissue_pvalues_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_tgfm_log_pvalue_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_pvalue_bar_plot, file=output_file, width=7.2, height=3.7, units="in")
	tmp_tissues <-get_sig_tissues_tgfm_tissue_pvalues_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	valid_tissues <- c(valid_tissues,tmp_tissues)
}
valid_tissues = sort(unique(valid_tissues))
}

##########################################################
# Heatmap showing expected number of causal genes in each tissue-trait pair
##########################################################
if (FALSE) {
pip_threshs <- c(0.0,.01, .1, .3, .5, .7)
pip_threshs <- c(.3, .5, .7)

for (pip_iter in 1:length(pip_threshs)) {
	pip_thresh <- pip_threshs[pip_iter]
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	heatmap <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits(trait_names, trait_names_readable, method_version, tgfm_results_dir, tissue_names, pip_thresh, valid_tissues)
	output_file <- paste0(visualize_tgfm_dir, "expected_num_causal_genes_", pip_thresh, "_", method_version,"_heatmap.pdf")
	ggsave(heatmap, file=output_file, width=7.2, height=8.0, units="in")
}
}

##########################################################
# Swarm-plot showing gene-tissue PIPs colored by tissue group for each trait
##########################################################
if (FALSE) {
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
beeswarm_plot <- make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait(trait_names, trait_names_readable, method_version, tgfm_results_dir)
output_file <- paste0(visualize_tgfm_dir, "beeswarm_gene_tissue_pip_colored_by_tissue_group_", method_version,".pdf")
ggsave(beeswarm_plot, file=output_file, width=15.2, height=4.5, units="in")
}





##########################################################
# Bar plot showing expected number of causal gene-tissue pairs in each tissue
##########################################################
if (FALSE) {
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]
	print(trait_name)

	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_gene_tissue_pairs_", trait_name_readable, "_", method_version,".pdf")
	ggsave(tissue_bar_plot, file=output_file, width=7.2, height=3.7, units="in")
}
}


##########################################################
# Barplot with standard errors showing fraction of mediated components across traits
##########################################################
if (FALSE) {
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_average_expression_mediated_probability_se_barplot.pdf")
med_prob_se_barplot <- make_mediated_prob_se_barplot(trait_names, trait_names_readable, method_version, tgfm_results_dir)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")


##########################################################
# Barplot with standard errors showing fraction of high pip genetic elements from gene expression
##########################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_average_fraction_of_high_pip_genetic_elements_from_gene_expression_se_barplot.pdf")
med_prob_se_barplot <- make_fraction_of_genetic_elements_from_gene_expression_se_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=4.2, units="in")

##########################################################
# Barplot with standard errors showing expected fraction of high pip genetic elements from gene expression
##########################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_expected_fraction_of_genetic_elements_from_gene_expression_se_barplot.pdf")
med_prob_se_barplot <- make_expected_fraction_of_genetic_elements_from_gene_expression_se_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=4.2, units="in")

}

##########################################################
# Barplot with standard errors showing number high pip genetic elements (stratefied by variants vs gene-tissue pairs)
##########################################################
if (FALSE) {
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
pip_thresholds <- c("0.25", "0.5", "0.75")
for (pip_iter in 1:length(pip_thresholds)) {
	pip_threshold = pip_thresholds[pip_iter]
	output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_", pip_threshold, "_elements_se_barplot.pdf")
	med_prob_se_barplot <- make_number_of_high_pip_elements_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits)
	ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")
}
}

##########################################################
# Barplot showing number high pip genetic gene-tissue pairs
##########################################################
if (FALSE) {
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
pip_thresholds <- c("0.25", "0.5", "0.75", "0.9")
for (pip_iter in 1:length(pip_thresholds)) {
	pip_threshold = pip_thresholds[pip_iter]
	output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_", pip_threshold, "_gene_tissue_pairs_barplot.pdf")
	med_prob_se_barplot <- make_number_of_high_pip_gene_tissue_pairs_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits)
	ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")
}
}

if (FALSE) {
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_barplot_cross_2_thresholds.pdf")
pip_thresholds <- c("0.5", "0.9")
med_prob_se_barplot <- make_number_of_high_pip_cross_2_threshold_gene_tissue_pairs_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_thresholds, independent_traits)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")

method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_barplot_cross_3_thresholds.pdf")
pip_thresholds <- c("0.5", "0.75", "0.9")
med_prob_se_barplot <- make_number_of_high_pip_cross_3_threshold_gene_tissue_pairs_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_thresholds, independent_traits)
ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")
}

##########################################################
# Make heatmap showing expected number of causal genes in each tissue-trait pair for selected traits
##########################################################
if (FALSE) {
pip_thresh <- "0.5"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
selected_traits <- c("mental_NEUROTICISM", "disease_AID_ALL", "disease_HYPERTENSION_DIAGNOSED", "disease_ALLERGY_ECZEMA_DIAGNOSED", "biochemistry_VitaminD", "body_HEIGHTz", "body_BMIz", "biochemistry_LDLdirect", "disease_HYPOTHYROIDISM_SELF_REP", "blood_MONOCYTE_COUNT", "biochemistry_HbA1c")
heatmap <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_thresh, selected_traits)
output_file <- paste0(visualize_tgfm_dir, "expected_num_causal_genes_", pip_thresh, "_", method_version,"_selected_traits_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=5.0, units="in")

pip_thresh <- "0.7"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
selected_traits <- c("mental_NEUROTICISM", "disease_AID_ALL", "disease_HYPERTENSION_DIAGNOSED", "disease_ALLERGY_ECZEMA_DIAGNOSED", "biochemistry_VitaminD", "body_HEIGHTz", "body_BMIz", "biochemistry_LDLdirect", "disease_HYPOTHYROIDISM_SELF_REP", "blood_MONOCYTE_COUNT", "biochemistry_HbA1c")
heatmap <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_thresh, selected_traits)
output_file <- paste0(visualize_tgfm_dir, "expected_num_causal_genes_", pip_thresh, "_", method_version,"_selected_traits_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=5.0, units="in")
}

##########################################################
# Make heatmap-barplot showing expected number of causal genes
##########################################################

method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file_png <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_heatmap_barplot.png")
output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_heatmap_barplot.pdf")
heatmap_barplot <- make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
ggsave(heatmap_barplot, file=output_file_png, width=7.2, height=3.7, units="in", dpi=400)
ggsave(heatmap_barplot, file=output_file_pdf, width=7.2, height=3.7, units="in", dpi=400)

method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file_png <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2.png")
output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2.pdf")
heatmap_barplot <- make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
ggsave(heatmap_barplot, file=output_file_png, width=7.2, height=3.7, units="in", dpi=400)
ggsave(heatmap_barplot, file=output_file_pdf, width=7.2, height=3.7, units="in", dpi=400)




if (FALSE) {
##########################################################
# Make Figure 3
##########################################################
pip_threshold="0.5"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"

#figure_3a <- plot_grid(NULL, make_number_of_high_pip_elements_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits), ncol=2, rel_widths=c(.03,1))
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
pip_thresholds <- c("0.5", "0.75", "0.9")
figure_3a <- plot_grid(NULL,make_number_of_high_pip_cross_3_threshold_gene_tissue_pairs_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_thresholds, independent_traits), ncol=2, rel_widths=c(.03,1))


selected_traits <- c("mental_NEUROTICISM", "disease_AID_ALL", "disease_HYPERTENSION_DIAGNOSED", "disease_ALLERGY_ECZEMA_DIAGNOSED", "biochemistry_VitaminD", "body_HEIGHTz", "body_BMIz", "biochemistry_LDLdirect", "disease_HYPOTHYROIDISM_SELF_REP", "blood_MONOCYTE_COUNT", "biochemistry_HbA1c")
figure_3b <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_threshold, selected_traits)

figure3 <- plot_grid(figure_3a, figure_3b, ncol=1, rel_heights=c(.65, 1.0), labels=c("a","b"))

output_file <- paste0(visualize_tgfm_dir, "figure_3.pdf")
ggsave(figure3, file=output_file, width=7.2, height=7.5, units="in")
}











print("DONE")




if (FALSE) {
##########################################################
# Scatter with standard errors showing fraction of mediated components across traits
##########################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_average_expression_mediated_vs_expression_mediated_h2_scatterplot.pdf")
med_prob_se_scatter <- make_mediated_prob_se_scatterplot(trait_names, trait_names_readable, method_version, tgfm_results_dir, tgfm_sldsc_results_dir)
ggsave(med_prob_se_scatter, file=output_file, width=7.2, height=4.7, units="in")
}


if (FALSE) {
#########################
# Make figure 2

# Panel 2a
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
panel_2a <- plot_grid(NULL, make_mediated_prob_se_barplot(trait_names, trait_names_readable, method_version, tgfm_results_dir), ncol=2, rel_widths=c(.015,1))


# Panel 2b
trait_name = "biochemistry_LDLdirect"
trait_name_readable = "LDL"
panel_2b <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait_for_figure(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)


# Panel 2c
trait_name = "disease_HYPERTENSION_DIAGNOSED"
trait_name_readable = "Hypertension"
panel_2c <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait_for_figure(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)


# Make 2bc panel
legender = get_legend(panel_2b)
panel_2bc <- plot_grid(plot_grid(panel_2b + theme(legend.position="none"), panel_2c + theme(legend.position="none"), ncol=2, labels=c("b", "c")), legender, ncol=1, rel_heights=c(1,.1))


pip_thresh = .3
print("START")
panel_2d <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits_for_figure2(trait_names, trait_names_readable, method_version, tgfm_results_dir, tissue_names, pip_thresh, valid_tissues)



figure2 <- plot_grid(panel_2a, panel_2bc, panel_2d, ncol=1, labels=c("a", "", "d"), rel_heights=c(.28,.24,.65))

output_file <- paste0(visualize_tgfm_dir, "figure2.pdf")
ggsave(figure2, file=output_file, width=7.2, height=11.5, units="in")


# Panel 2a
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
panel_2a <- plot_grid(NULL, make_mediated_prob_se_barplot(trait_names, trait_names_readable, method_version, tgfm_results_dir), ncol=2, rel_widths=c(.015,1))


# Panel 2b
trait_name = "biochemistry_LDLdirect"
trait_name_readable = "LDL"
panel_2b <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait_for_figure(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)


# Panel 2c
trait_name = "disease_HYPERTENSION_DIAGNOSED"
trait_name_readable = "Hypertension"
panel_2c <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait_for_figure(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)


# Make 2bc panel
legender = get_legend(panel_2b)
panel_2bc <- plot_grid(plot_grid(panel_2b + theme(legend.position="none"), panel_2c + theme(legend.position="none"), ncol=2, labels=c("b", "c")), legender, ncol=1, rel_heights=c(1,.1))


pip_thresh = .5
print("START")
panel_2d <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits_for_figure2(trait_names, trait_names_readable, method_version, tgfm_results_dir, tissue_names, pip_thresh, valid_tissues)



figure2 <- plot_grid(panel_2a, panel_2bc, panel_2d, ncol=1, labels=c("a", "", "d"), rel_heights=c(.25,.22,.66))

output_file <- paste0(visualize_tgfm_dir, "figure2_v2.pdf")
ggsave(figure2, file=output_file, width=7.2, height=10.5, units="in")



# Panel 2a
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
panel_2a <- plot_grid(NULL, make_mediated_prob_se_barplot(trait_names, trait_names_readable, method_version, tgfm_results_dir), ncol=2, rel_widths=c(.015,1))


# Panel 2b
trait_name = "biochemistry_LDLdirect"
trait_name_readable = "LDL"
panel_2b <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait_for_figure(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)


# Panel 2c
trait_name = "disease_HYPERTENSION_DIAGNOSED"
trait_name_readable = "Hypertension"
panel_2c <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait_for_figure(trait_name, trait_name_readable, method_version, tgfm_results_dir, tissue_names)


# Make 2bc panel
legender = get_legend(panel_2b)
panel_2bc <- plot_grid(plot_grid(panel_2b + theme(legend.position="none"), panel_2c + theme(legend.position="none"), ncol=2, labels=c("b", "c")), legender, ncol=1, rel_heights=c(1,.1))


pip_thresh = .7
print("START")
panel_2d <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits_for_figure2(trait_names, trait_names_readable, method_version, tgfm_results_dir, tissue_names, pip_thresh, valid_tissues)



figure2 <- plot_grid(panel_2a, panel_2bc, panel_2d, ncol=1, labels=c("a", "", "d"), rel_heights=c(.25,.22,.66))

output_file <- paste0(visualize_tgfm_dir, "figure2_v3.pdf")
ggsave(figure2, file=output_file, width=7.2, height=10.5, units="in")
}






if (FALSE) {
##########################################################
# Component level tissue plot
##########################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_average_expression_mediated_probability_per_tissue_se_barplot.pdf")
med_prob_se_barplot <- make_tissue_mediated_prob_se_barplot(trait_names, trait_names_readable, method_version, tgfm_results_dir)
#ggsave(med_prob_se_barplot, file=output_file, width=7.2, height=3.7, units="in")
}

if (FALSE) {

##########################################################
# Bar plot showing expected number of causal tissue categories
##########################################################
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	#tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_tissue_categories_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_results_dir)
	#output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_tissue_categories_", trait_name_readable, "_", method_version,".pdf")
	#ggsave(tissue_bar_plot, file=output_file, width=7.2, height=4.5, units="in")

}

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


