args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(reshape2)
library(ggbeeswarm)
library(RColorBrewer)
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

make_number_of_high_pip_genes_heatmap_barplot_v2_transpose <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, gene_type="component_gene",preordered=FALSE, ordered_traits=NULL) {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_tgfm_n_causal_genes_cross_pip_threshold_sqrt_plot_input.txt")
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

	indices = order(mid_point_vec, total_hits_vec)

	if (preordered == FALSE) {
		df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)
	} else {
		df$trait = factor(df$trait, levels=ordered_traits)
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=ordered_traits), midpoint=mid_point_vec)
	}

   p <- ggplot() +
  	geom_bar(data=df, aes(y = trait,x = num_elements,fill = PIP, color=PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	scale_x_continuous(breaks=c(0.0,sqrt(25), sqrt(100), sqrt(250), sqrt(500), sqrt(1000)), labels=c("0", "25", "100", "250", "  500", "1000")) +
  	#theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_distiller(direction=1, palette = "Greens", limits=c(.2,1.0)) +
  	scale_color_distiller(direction=1, palette = "Greens", limits=c(.2,1.0)) +	
  	theme(plot.title = element_text(hjust = 0.5)) +
  	 geom_errorbar(data  = df2, aes(y=trait,x=midpoint, xmax=midpoint, xmin=midpoint)) +
  	 labs(y="", x="No. fine-mapped genes")
 	return(p)
}

make_number_of_high_pip_genes_heatmap_barplot_v2 <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, preordered=FALSE, ordered_traits=NULL) {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genes_cross_pip_threshold_sqrt_plot_input.txt")
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

	if (preordered == FALSE) {
		df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)
	} else {
		df$trait = factor(df$trait, levels=ordered_traits)
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=ordered_traits), midpoint=mid_point_vec)
	}

   p <- ggplot() +
  	geom_bar(data=df, aes(x = trait,y = num_elements,fill = PIP, color=PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	scale_y_continuous(breaks=c(0.0,sqrt(25), sqrt(100), sqrt(250), sqrt(500), sqrt(1000)), labels=c("0", "25", "100", "250", "  500", "1000")) +
  	theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_distiller(direction=1, palette = "Greens", limits=c(.2,1.0)) +
  	scale_color_distiller(direction=1, palette = "Greens", limits=c(.2,1.0)) +	
  	theme(plot.title = element_text(hjust = 0.5)) +
  	 geom_errorbar(data  = df2, aes(x=trait,y=midpoint, ymax=midpoint, ymin=midpoint)) +
  	 labs(x="", y="No. fine-mapped\nGenes")
 	return(p)
}

make_number_of_high_pip_variants_heatmap_barplot_v2_transpose <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, gene_type="component_gene", preordered=FALSE, ordered_traits=NULL) {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_", gene_type,"_", method_version, "_tgfm_n_causal_variants_cross_pip_threshold_sqrt_plot_input.txt")
			tmp_data <- read.table(summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="variant",]

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


	indices = order(mid_point_vec,total_hits_vec)

	if (preordered == FALSE) {
		df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)
	} else {
		df$trait = factor(df$trait, levels=ordered_traits)
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=ordered_traits), midpoint=mid_point_vec)
	}

   p <- ggplot() +
  	geom_bar(data=df, aes(y = trait,x = num_elements,fill = PIP, color=PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	scale_x_continuous(breaks=c(0.0, sqrt(100), sqrt(500), sqrt(1000), sqrt(2000)), labels=c(0, 100, 500, 1000, 2000)) +
  	#theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_distiller(direction=1, palette = "Blues", limits=c(.2,1.0)) +
  	scale_color_distiller(direction=1, palette = "Blues", limits=c(.2,1.0)) +	
  	theme(plot.title = element_text(hjust = 0.5)) +
  	 geom_errorbar(data  = df2, aes(y=trait,x=midpoint, xmax=midpoint, xmin=midpoint)) +
  	 labs(y="", x="No. fine-mapped variants")
 	return(p)
}



make_number_of_high_pip_variants_heatmap_barplot_v2 <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, preordered=FALSE, ordered_traits=NULL) {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_variants_cross_pip_threshold_sqrt_plot_input.txt")
			tmp_data <- read.table(summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="variant",]

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

	if (preordered == FALSE) {
		df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)
	} else {
		df$trait = factor(df$trait, levels=ordered_traits)
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=ordered_traits), midpoint=mid_point_vec)
	}

   p <- ggplot() +
  	geom_bar(data=df, aes(x = trait,y = num_elements,fill = PIP, color=PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	scale_y_continuous(breaks=c(0.0, sqrt(100), sqrt(500), sqrt(1000), sqrt(2000)), labels=c(0, 100, 500, 1000, 2000)) +
  	theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_distiller(direction=1, palette = "Blues", limits=c(.2,1.0)) +
  	scale_color_distiller(direction=1, palette = "Blues", limits=c(.2,1.0)) +	
  	theme(plot.title = element_text(hjust = 0.5)) +
  	 geom_errorbar(data  = df2, aes(x=trait,y=midpoint, ymax=midpoint, ymin=midpoint)) +
  	 labs(x="", y="No. fine-mapped\nVariants")
 	return(p)
}

make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_trait_order_according_to_gene_tissue_pairs<- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, gene_type="component_gene") {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_", gene_type,"_", method_version, "_tgfm_n_causal_gene_tissue_pairs_cross_pip_threshold_sqrt_plot_input.txt")
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

	indices = order(-mid_point_vec)
	df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])

	return(as.character(trait_names_vec2)[indices])
}

make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2_transpose <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, preordered=FALSE, ordered_traits=NULL, gene_type="component_gene") {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_tgfm_n_causal_gene_tissue_pairs_cross_pip_threshold_sqrt_plot_input.txt")
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

	indices = order(mid_point_vec, total_hits_vec)

	if (preordered == FALSE) {
		df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)
	} else {
		df$trait = factor(df$trait, levels=ordered_traits)
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=ordered_traits), midpoint=mid_point_vec)		
	}

   p <- ggplot() +
  	geom_bar(data=df, aes(y = trait,x = num_elements,fill = PIP, color=PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	#scale_y_continuous(breaks=c(0.0,sqrt(5), sqrt(20), sqrt(50), sqrt(100), sqrt(200), sqrt(400), sqrt(600)), labels=c(0, 5, 20, 50, 100, 200, 400,600)) +
  	scale_x_continuous(breaks=c(0.0,sqrt(10), sqrt(50), sqrt(100), sqrt(200)), labels=c("0", "10", "50", "100", "  200")) +
  	#theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_distiller(direction=1, palette = "Reds", limits=c(.2,1.0)) +
  	scale_color_distiller(direction=1, palette = "Reds", limits=c(.2,1.0)) +	
  	theme(plot.title = element_text(hjust = 0.5)) +
  	 geom_errorbar(data  = df2, aes(y=trait,x=midpoint, xmax=midpoint, xmin=midpoint)) +
  	 labs(y="", x="No. fine-mapped gene-tissue pairs")
 	return(p)
}


make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2 <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, preordered=FALSE, ordered_traits=NULL) {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_gene_tissue_pairs_cross_pip_threshold_sqrt_plot_input.txt")
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

	if (preordered == FALSE) {
		df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]), midpoint=mid_point_vec)
	} else {
		df$trait = factor(df$trait, levels=ordered_traits)
		df2 <- data.frame(trait=factor(trait_names_vec2, levels=ordered_traits), midpoint=mid_point_vec)		
	}

   p <- ggplot() +
  	geom_bar(data=df, aes(x = trait,y = num_elements,fill = PIP, color=PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	#scale_y_continuous(breaks=c(0.0,sqrt(5), sqrt(20), sqrt(50), sqrt(100), sqrt(200), sqrt(400), sqrt(600)), labels=c(0, 5, 20, 50, 100, 200, 400,600)) +
  	scale_y_continuous(breaks=c(0.0,sqrt(10), sqrt(50), sqrt(100), sqrt(200)), labels=c("0", "10", "50", "100", "  200")) +
  	theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_distiller(direction=1, palette = "Reds", limits=c(.2,1.0)) +
  	scale_color_distiller(direction=1, palette = "Reds", limits=c(.2,1.0)) +	
  	theme(plot.title = element_text(hjust = 0.5)) +
  	 geom_errorbar(data  = df2, aes(x=trait,y=midpoint, ymax=midpoint, ymin=midpoint)) +
  	 labs(x="", y="No. fine-mapped\nGene-Tissue pairs")
 	return(p)
}

make_number_of_high_pip_gene_tissue_pairs_categorical_heatmap_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits) {
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_5_pip_bins_sqrt_plot_input.txt")
			tmp_data <- read.table(summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="gene",]
			
			number_elements_vec <- c(number_elements_vec, tmp_data$n_elements)
			pip_threshold_vec <- c(pip_threshold_vec, as.character(tmp_data$PIP_threshold))
			trait_names_vec <- c(trait_names_vec, rep(trait_name_readable, length(tmp_data$PIP_threshold)))

			total_hits_vec <- c(total_hits_vec, sum(tmp_data$n_elements))
			trait_names_vec2 <- c(trait_names_vec2, trait_name_readable)
		}
	}
	df <- data.frame(trait=trait_names_vec, PIP=factor(pip_threshold_vec), num_elements=number_elements_vec)
	indices = order(-total_hits_vec)
	df$trait = factor(df$trait, levels=as.character(trait_names_vec2)[indices])


	#df2 <- data.frame(trait=factor(trait_names_vec2, levels=as.character(trait_names_vec2)[indices]))
	blues <- brewer.pal(5, "Blues")
	blue_range<-colorRampPalette(blues)


   p <- ggplot() +
  	geom_bar(data=df, aes(x = trait,y = num_elements,fill = PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	scale_y_continuous(breaks=c(0.0,sqrt(5), sqrt(20), sqrt(50), sqrt(100), sqrt(200), sqrt(400), sqrt(600)), labels=c(0, 5, 20, 50, 100, 200, 400,600)) +
  	#theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
  	theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_manual(values = blue_range(5)) +
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

			summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_cross_pip_threshold_sqrt_plot_input.txt")
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


	tot_var_count = 0
	tot_gene_count = 0
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {

			summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_expected_n_causal_genetic_elements.txt")
			tmp <- read.table(summary_file, header=TRUE)


			variant_count = tmp$count[1]
			gene_count = tmp$count[2]
			total_count = gene_count + variant_count

			tot_var_count = tot_var_count + variant_count
			tot_gene_count = tot_gene_count + gene_count
			mean_mediated_prob <- gene_count/(total_count)
		
			se_mediated_prob = sqrt((mean_mediated_prob*(1.0-mean_mediated_prob))/(total_count))
			trait_names_vec <- c(trait_names_vec, trait_name_readable)
			fraction_mediated_vec <- c(fraction_mediated_vec, mean_mediated_prob)
			fraction_mediated_se_vec <- c(fraction_mediated_se_vec, se_mediated_prob)
		}
	}


	df <- data.frame(trait=trait_names_vec, average_mediated_probability=fraction_mediated_vec, average_mediated_probability_se=fraction_mediated_se_vec)
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	indices = order(-df$average_mediated_probability)
	df$trait = factor(df$trait, levels=as.character(df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=average_mediated_probability)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		#scale_fill_manual(values=c("plum2", "orchid2", "magenta2")) +
    		theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
    		labs(y="Gene-mediated\nPIP fraction", x="") +
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

	gene_tot_count = 0
	total_total_count = 0


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

			if (pip_threshold == "0.5") {
				gene_tot_count = gene_tot_count + gene_count
				total_total_count = total_total_count + gene_count + variant_count
			}
		}


		}
	}

	df <- data.frame(trait=trait_names_vec, average_mediated_probability=fraction_mediated_vec, average_mediated_probability_se=fraction_mediated_se_vec, pip=factor(pip_threshold_vec,levels=pip_thresholds))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$pip) == "0.25",]

	indices = order(-tmp_df$average_mediated_probability)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	p <- ggplot(df, aes(x=trait, y=average_mediated_probability, fill=pip)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		scale_fill_manual(values=c("plum2", "orchid2", "magenta2")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="Proportion of\n fine-mapped\ngenetic elements", x="") +
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
bh_fdr_correction <- function(p_values, alpha=0.05) {
  # Get the number of p-values
  m <- length(p_values)
  
  # Create an index vector for sorting
  index <- order(p_values)
  
  # Sort the p-values in ascending order
  sorted_p_values <- p_values[index]
  
  # Calculate the Benjamini-Hochberg critical values
  critical_values <- (1:m) * (alpha / m)
  
  # Create a vector of corrected p-values
  adjusted_p_values <- rep(NA, m)

  if (sum(sorted_p_values <= critical_values) > 0) {

 	 # Identify the largest k such that p_k <= critical_values[k]
 	 k <- max(which(sorted_p_values <= critical_values))

 	 adjusted_p_values[index[1:k]] <- sorted_p_values[1:k]
	}
  
  return(adjusted_p_values)
}



generate_file_containing_bonf_significance_of_each_trait_tissue_pair_based_on_iterative_prior <- function(trait_names, iterative_tgfm_prior_dir, method_version,trait_tissue_prior_significance_file, gene_type) {
	trait_vec <- c()
	tissue_vec <- c()
	sig_vec <- c()
	nom_vec <- c()
	fdr_sig_vec <- c()
	mean_prob_vec <- c()
	se_prob_vec <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name = trait_names[trait_iter]
		iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
	
		df <- read.table(iterative_prior_file, header=TRUE)
		# Skip variants
		nrow = dim(df)[1]

		df <- df[2:nrow,]
		df$tissue = df$element_name

		n_tissues = length(unique(df$tissue))

		trait_nom_vec <- c()
		trait_mean_prob_vec <- c()
		for (tissue_iter in 1:length(df$tissue)) {
			tissue_name <- as.character(df$tissue[tissue_iter])
			prob_string = as.character(df$prior_distribution[tissue_iter])
			prob_vec = as.numeric(strsplit(prob_string,";")[[1]])

			sdev = sd(prob_vec)
			if (sdev == 0.0) {
				zz = 0
			} else {
				zz = mean(prob_vec)/(sd(prob_vec))
			}

			pvalue = pnorm(q=abs(zz), lower.tail=FALSE)
			bonf_pvalue = pvalue*n_tissues

			trait_vec <- c(trait_vec, trait_name)
			tissue_vec <- c(tissue_vec, tissue_name)
			sig_vec <- c(sig_vec, bonf_pvalue)
			nom_vec <- c(nom_vec, pvalue)
			trait_nom_vec <- c(trait_nom_vec, pvalue)
			mean_prob_vec <- c(mean_prob_vec, mean(prob_vec))
			se_prob_vec <- c(se_prob_vec, sd(prob_vec))
			trait_mean_prob_vec <- c(trait_mean_prob_vec,mean(prob_vec) )
		}
		bh_corrected_pvalues_05 = bh_fdr_correction(trait_nom_vec, alpha=0.05)
		bh_corrected_pvalues_2 = bh_fdr_correction(trait_nom_vec, alpha=0.2)


		for (tissue_iter in 1:length(bh_corrected_pvalues_05)) {
			if (is.na(bh_corrected_pvalues_05[tissue_iter]) == FALSE) {
				if (trait_mean_prob_vec[tissue_iter] > 1e-14) {
					fdr_sig_vec <- c(fdr_sig_vec, "**")
				} else {
					fdr_sig_vec <- c(fdr_sig_vec, "null")
				}
			} else if (is.na(bh_corrected_pvalues_2[tissue_iter]) == FALSE) {
				if (trait_mean_prob_vec[tissue_iter] > 1e-14) {
					fdr_sig_vec <- c(fdr_sig_vec, "*")
				} else {
					fdr_sig_vec <- c(fdr_sig_vec, "null")
				}
			} else {
				fdr_sig_vec <- c(fdr_sig_vec, "null")
			}
		}


	}

	df <- data.frame(tissue=tissue_vec, trait=trait_vec, mean_prob=mean_prob_vec, se_prob=se_prob_vec, pvalue=nom_vec, fdr_significance=fdr_sig_vec)
	write.table(df, file=trait_tissue_prior_significance_file, quote=FALSE, sep="\t", row.names = FALSE)
}



generate_file_containing_bonf_significance_of_each_trait_tissue_pair_based_on_iterative_prior_alt <- function(trait_names, iterative_tgfm_prior_dir, method_version,trait_tissue_prior_significance_file, gene_type, fdr_thresholds) {
	trait_vec <- c()
	tissue_vec <- c()
	sig_vec <- c()
	nom_vec <- c()
	fdr_sig_vec <- c()
	mean_prob_vec <- c()
	se_prob_vec <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name = trait_names[trait_iter]
		iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
	
		df <- read.table(iterative_prior_file, header=TRUE)
		# Skip variants
		nrow = dim(df)[1]

		df <- df[2:nrow,]
		df$tissue = df$element_name

		n_tissues = length(unique(df$tissue))

		trait_nom_vec <- c()
		trait_mean_prob_vec <- c()
		for (tissue_iter in 1:length(df$tissue)) {
			tissue_name <- as.character(df$tissue[tissue_iter])
			prob_string = as.character(df$prior_distribution[tissue_iter])
			prob_vec = as.numeric(strsplit(prob_string,";")[[1]])

			sdev = sd(prob_vec)
			if (sdev == 0.0) {
				zz = 0
			} else {
				zz = mean(prob_vec)/(sd(prob_vec))
			}

			pvalue = pnorm(q=abs(zz), lower.tail=FALSE)
			bonf_pvalue = pvalue*n_tissues

			trait_vec <- c(trait_vec, trait_name)
			tissue_vec <- c(tissue_vec, tissue_name)
			sig_vec <- c(sig_vec, bonf_pvalue)
			nom_vec <- c(nom_vec, pvalue)
			trait_nom_vec <- c(trait_nom_vec, pvalue)
			mean_prob_vec <- c(mean_prob_vec, mean(prob_vec))
			se_prob_vec <- c(se_prob_vec, sd(prob_vec))
			trait_mean_prob_vec <- c(trait_mean_prob_vec,mean(prob_vec) )
		}

		ordered_fdr_thresholds = rev(sort(fdr_thresholds))

		ignored_indices = trait_mean_prob_vec <= 1e-14

		trait_fdr_sig_vec <- rep("null", length(trait_mean_prob_vec))

		for (fdr_thresh_iter in 1:length(ordered_fdr_thresholds)) {
			fdr_threshold <- ordered_fdr_thresholds[fdr_thresh_iter]
			bh_corrected_pvalues = bh_fdr_correction(trait_nom_vec, alpha=fdr_threshold)
			print(bh_corrected_pvalues)
		}


		bh_corrected_pvalues_05 = bh_fdr_correction(trait_nom_vec, alpha=0.05)
		bh_corrected_pvalues_2 = bh_fdr_correction(trait_nom_vec, alpha=0.2)


		for (tissue_iter in 1:length(bh_corrected_pvalues_05)) {
			if (is.na(bh_corrected_pvalues_05[tissue_iter]) == FALSE) {
				if (trait_mean_prob_vec[tissue_iter] > 1e-14) {
					fdr_sig_vec <- c(fdr_sig_vec, "**")
				} else {
					fdr_sig_vec <- c(fdr_sig_vec, "null")
				}
			} else if (is.na(bh_corrected_pvalues_2[tissue_iter]) == FALSE) {
				if (trait_mean_prob_vec[tissue_iter] > 1e-14) {
					fdr_sig_vec <- c(fdr_sig_vec, "*")
				} else {
					fdr_sig_vec <- c(fdr_sig_vec, "null")
				}
			} else {
				fdr_sig_vec <- c(fdr_sig_vec, "null")
			}
		}


	}

	df <- data.frame(tissue=tissue_vec, trait=trait_vec, mean_prob=mean_prob_vec, se_prob=se_prob_vec, pvalue=nom_vec, fdr_significance=fdr_sig_vec)
	write.table(df, file=trait_tissue_prior_significance_file, quote=FALSE, sep="\t", row.names = FALSE)
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

		zz = mean(prob_vec)/(sd(prob_vec))

		pvalue = pnorm(q=abs(zz), lower.tail=FALSE)
		#pvalue = sum(prob_vec < threshold)/length(prob_vec)

		if (is.na(pvalue)) {
			pvalue = 1.0
		}

		if (pvalue < (.05/38)) {
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

make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh, selected_traits, trait_tissue_prior_significance_file, significance_bool=FALSE, only_fdr_05=FALSE) {
	df_prior = read.table(trait_tissue_prior_significance_file, header=TRUE)
	df_prior$fdr_significance = as.character(df_prior$fdr_significance)

	tissue_vec <- c()
	count_vec <- c()
	trait_vec <- c()
	sig_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		if (trait_name %in% selected_traits) {
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name_readable == "Bone mineral density") {
			trait_name_readable = "Bone density"
		}

		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")
		tmper <- c()
		tmp_df <- trait_df[trait_df$PIP > as.numeric(pip_thresh),]
		df_prior_trait_sub = df_prior[as.character(df_prior$trait)==trait_name,]
		if (sum(tmp_df$PIP) > 2.0) {
			#print(head(tmp_df))
			for (tissue_iter in 1:length(ordered_tissues)) {
				tissue_name <- ordered_tissues[tissue_iter]
				#pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				pip_sum = length(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				#count_vec <- c(count_vec, pip_sum)
				trait_vec <- c(trait_vec, trait_name_readable)
				tmper <- c(tmper, pip_sum)

				# Extract bonferonni significane from prior
				if (FALSE) {
				bonf_pvalue = df_prior_trait_sub$bonferonni_pvalue[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				if (bonf_pvalue <= .05) {
					sig_vec <- c(sig_vec, "**")
				} else if (bonf_pvalue <= .2) {
					sig_vec <- c(sig_vec, "*")
				} else {
					sig_vec <- c(sig_vec, "")
				}
				}
				sig_value <- df_prior_trait_sub$fdr_significance[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				if (sig_value == "null") {
					sig_vec <- c(sig_vec, "")
				} else {
					if (only_fdr_05) {
						if (sig_value == "**") {
							sig_vec <- c(sig_vec, "*")
						} else {
							sig_vec <- c(sig_vec, "")
						}
					} else {
						sig_vec <- c(sig_vec, sig_value)
						}
				}

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
	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues), trait=factor(trait_vec, levels=new_trait_names), significance=sig_vec)

	for (tissue_iter in 1:length(ordered_tissues)) {
		tissue_name = ordered_tissues[tissue_iter]
		indices = as.character(df$tissue) == tissue_name
		if (max(df$expected_causal_genes[indices]) > .2) {
			valid_tissues <- c(valid_tissues, tissue_name)
		}

	}
	valid_tissues <- c(valid_tissues, "Lung")


	df = df[as.character(df$tissue) %in% valid_tissues,]

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Spleen="spleen", Pituitary="pituitary", Liver="liver",Lung="lung", Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="adipose visc.", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="skin (sun)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain Basal Ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord", Whole_Blood="whole blood", Colon_Sigmoid="Colon Sigmoid", Artery_Tibial="artery tibial", Artery_Aorta="artery aorta", Brain_Cerebellum="cerebellum", Brain_BasalGanglia="brain basal ganglia", Esophagus_Mucosa="esoph. mucosa", Brain_Cortex="brain cortex", Artery_Heart="artery heart", Adrenal_Gland="adrenal gland")
	df$tissue = factor(df$tissue)
	df$value = df$expected_causal_genes




	matrix <- dcast(df, tissue ~ trait)
	tissue_names <- matrix$tissue
	matrix <- as.matrix(matrix[,2:dim(matrix)[2]])

	ord <- hclust( dist(matrix, method = "euclidean"), method = "ward.D" )$order
	df$tissue <- factor(df$tissue, levels=as.character(tissue_names)[ord])
	custom_tissue_names <- c("spleen", "lymphocytes","whole blood", "skin (sun)", "esoph. mucosa", "brain basal ganglia","cerebellum", "pituitary", "adipose visc.", "liver", "artery aorta","artery tibial", "adrenal gland", "lung", "fibroblast", "thyroid", "Brain_Limbic")
	df$tissue <- factor(df$tissue, levels=custom_tissue_names)
	
	#ord2 <- hclust( dist(t(matrix), method = "euclidean"), method = "ward.D" )$order
	#df$trait <- factor(df$trait, levels=as.character(new_trait_names)[ord2])
	custom_orderd_traits = rev(c("All autoimmune", "Monocyte count", "Corp. hemoglobin", "Platelet volume", "Reticulocyte count", "Eczema", "Vitamin D", "Balding", "Hair color", "Menarche age", "Cholesterol", "Diastolic BP", "FVC", "FEV1:FVC", "Height", "Bone density", "Chronotype"))
	df$trait <- factor(df$trait, levels=custom_orderd_traits)


	red_color =brewer.pal(n = 9, name = "Reds")[7]

	#print(df)

	pp <- ggplot(df, aes(tissue, trait, fill= expected_causal_genes)) + 
  		geom_tile() +
  		figure_theme() +
  		#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low = "white", high = red_color) +
  		labs(x="",y="", fill="Fraction of fine-mapped gene-tissue pairs")

  	if (significance_bool == TRUE) {
  		pp = pp + geom_text(aes(label=significance),vjust=.8)
  	}

  	return(pp)	
}


make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits_old <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh, selected_traits) {
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


	red_color =brewer.pal(n = 9, name = "Reds")[7]

	print(df)

	pp <- ggplot(df, aes(tissue, trait, fill= expected_causal_genes)) + 
  		geom_tile() +
  		figure_theme() +
  		#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low = "white", high = red_color) +
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
		if (sum(trait_df$PIP > pip_thresh) > 0) {
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

get_heatmap_data_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh, trait_tissue_prior_significance_file) {
	df_prior = read.table(trait_tissue_prior_significance_file, header=TRUE)
	df_prior$fdr_significance = as.character(df_prior$fdr_significance)
	tissue_vec <- c()
	count_vec <- c()
	trait_vec <- c()
	sig_vec <- c()
	pvalue_vec <- c()
	bonf_pvalue_vec <- c()



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_full_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")

		df_prior_trait_sub = df_prior[as.character(df_prior$trait)==trait_name,]
		tmper <- c()
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		if (sum(trait_df$PIP > pip_thresh) > 0) {
			#print(head(tmp_df))
			for (tissue_iter in 1:length(ordered_tissues)) {
				tissue_name <- ordered_tissues[tissue_iter]
				pip_sum = length(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				#count_vec <- c(count_vec, pip_sum)
				trait_vec <- c(trait_vec, trait_name_readable)
				tmper <- c(tmper, pip_sum)
				# Extract bonferonni significane from prior
				#bonf_pvalue = df_prior_trait_sub$bonferonni_pvalue[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				pvalue = df_prior_trait_sub$pvalue[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				if (FALSE) {
				if (bonf_pvalue <= .05) {
					sig_vec <- c(sig_vec, "**")
				} else if (bonf_pvalue <= .2) {
					sig_vec <- c(sig_vec, "*")
				} else {
					sig_vec <- c(sig_vec, "")
				}
				}
				sig_value <- df_prior_trait_sub$fdr_significance[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				if (sig_value == "null") {
					sig_vec <- c(sig_vec, "null")
				} else {
					sig_vec <- c(sig_vec, sig_value)
				}


				pvalue_vec <- c(pvalue_vec, pvalue)
				#bonf_pvalue_vec <- c(bonf_pvalue_vec, bonf_pvalue)
			}
			tmper = tmper/sum(tmper)
			for (itera in 1:length(tmper)) {
				tissue_name <- ordered_tissues[itera]
				count_vec <- c(count_vec, tmper[itera])
			}
		}
	}

	new_trait_names = sort(unique(trait_vec))

	df <- data.frame(trait=factor(trait_vec, levels=new_trait_names), tissue=factor(tissue_vec,levels=ordered_tissues),proportion_fine_mapped_gene_tissue_pairs=count_vec, pvalue=pvalue_vec, fdr_significance=sig_vec)


  	return(df)


}

make_tissue_tissue_correlation_heatmap <- function(tissue_tissue_correlation_file) {
	df <- read.table(tissue_tissue_correlation_file, header=TRUE)

	df$value = df$average_correlation

	df$tissue_ii = str_replace_all(as.character(df$tissue_ii), "-", "_")
	df$tissue_ii <- recode(df$tissue_ii, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$tissue_jj = str_replace_all(as.character(df$tissue_jj), "-", "_")
	df$tissue_jj <- recode(df$tissue_jj, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")



	matrix <- dcast(df, tissue_ii ~ tissue_jj)
	tissue_names <- matrix$tissue_ii
	matrix <- as.matrix(matrix[,2:dim(matrix)[2]])

	ord <- hclust( dist(matrix, method = "euclidean"), method = "ward.D" )$order
	print(as.character(tissue_names))
	print(as.character(tissue_names)[ord])

	df$tissue_ii <- factor(df$tissue_ii, levels=as.character(tissue_names)[ord])
	df$tissue_jj <- factor(df$tissue_jj, levels=as.character(tissue_names)[ord])

	pp <- ggplot(df, aes(tissue_ii, tissue_jj, fill= average_correlation)) + 
  		geom_tile() +
  		figure_theme() +
  		labs(fill="Avg. correlation",x="", y="") +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low = "white", high = "red")
  	return(pp)
}




make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh, trait_tissue_prior_significance_file, gene_type) {
	df_prior = read.table(trait_tissue_prior_significance_file, header=TRUE)
	df_prior$fdr_significance = as.character(df_prior$fdr_significance)
	tissue_vec <- c()
	count_vec <- c()
	trait_vec <- c()
	sig_vec <- c()
	pvalue_vec <- c()



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_tgfm_per_gene_tissue_full_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")

		df_prior_trait_sub = df_prior[as.character(df_prior$trait)==trait_name,]
		tmper <- c()
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		if (sum(trait_df$PIP > pip_thresh) > 0) {
			#print(head(tmp_df))
			for (tissue_iter in 1:length(ordered_tissues)) {
				tissue_name <- ordered_tissues[tissue_iter]
				pip_sum = length(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				#count_vec <- c(count_vec, pip_sum)
				trait_vec <- c(trait_vec, trait_name_readable)
				tmper <- c(tmper, pip_sum)
				# Extract bonferonni significane from prior
				#bonf_pvalue = df_prior_trait_sub$bonferonni_pvalue[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				pvalue = df_prior_trait_sub$pvalue[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				if (FALSE) {
				if (bonf_pvalue <= .05) {
					sig_vec <- c(sig_vec, "**")
				} else if (bonf_pvalue <= .2) {
					sig_vec <- c(sig_vec, "*")
				} else {
					sig_vec <- c(sig_vec, "")
				}
				}
				sig_value <- df_prior_trait_sub$fdr_significance[as.character(df_prior_trait_sub$tissue)==tissue_name][1]
				if (sig_value == "null") {
					sig_vec <- c(sig_vec, "")
				} else {
					sig_vec <- c(sig_vec, sig_value)
				}
				pvalue_vec <- c(pvalue_vec, pvalue)
				#bonf_pvalue_vec <- c(bonf_pvalue_vec, bonf_pvalue)
			}
			tmper = tmper/sum(tmper)
			for (itera in 1:length(tmper)) {
				tissue_name <- ordered_tissues[itera]
				count_vec <- c(count_vec, tmper[itera])
			}
		}
	}

	new_trait_names = sort(unique(trait_vec))

	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues), trait=factor(trait_vec, levels=new_trait_names),significance=sig_vec, pvalue=pvalue_vec)
	for (tissue_iter in 1:length(ordered_tissues)) {
		tissue_name = ordered_tissues[tissue_iter]
		indices = as.character(df$tissue) == tissue_name
	}

	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")
	df$tissue = factor(df$tissue)
	df$value = df$expected_causal_genes


	#df = df[as.character(df$tissue) %in% valid_tissues,]




	matrix <- dcast(df, tissue ~ trait)
	tissue_names <- matrix$tissue
	matrix <- as.matrix(matrix[,2:dim(matrix)[2]])

	ord <- hclust( dist(matrix, method = "euclidean"), method = "ward.D" )$order
	df$tissue <- factor(df$tissue, levels=as.character(tissue_names)[ord])
	ord2 <- hclust( dist(t(matrix), method = "euclidean"), method = "ward.D" )$order

	df$trait <- factor(df$trait, levels=as.character(new_trait_names)[ord2])

	pp <- ggplot(df, aes(tissue, trait, fill= expected_causal_genes)) + 
  		geom_tile() +
  		geom_text(aes(label=significance)) +
  		figure_theme() +
  		labs(fill="Proportion of fine-mapped gene-tissue pairs",x="", y="") +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  		theme(legend.position="bottom") +
  		scale_fill_gradient(low = "white", high = "red")

  	return(pp)


}


make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_normalized_by_number_of_genes_per_tissue_cross_traits <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, ordered_tissues, pip_thresh) {
	tissue_vec <- c()
	count_vec <- c()
	trait_vec <- c()


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_per_gene_tissue_full_pip_summary.txt")
		trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")


		tmper <- c()
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		if (sum(trait_df$PIP > pip_thresh) > 0) {
			#print(head(tmp_df))
			for (tissue_iter in 1:length(ordered_tissues)) {
				tissue_name <- ordered_tissues[tissue_iter]
				total_genes = sum(as.character(trait_df$tissue_name)==tissue_name)
				pip_sum = sum(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
				tissue_vec <- c(tissue_vec, tissue_name)
				#count_vec <- c(count_vec, pip_sum)
				trait_vec <- c(trait_vec, trait_name_readable)
				tmper <- c(tmper, pip_sum/total_genes)
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


	#df = df[as.character(df$tissue) %in% valid_tissues,]



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
  		scale_fill_gradient(low = "white", high = "red")

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


make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait <- function(trait_name, trait_name_readable, method_version, tgfm_results_dir, ordered_tissues, gene_type) {
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

	pip_threshs <- c(.5)
	#pip_threshs <- c(.3, .5, .7)

	trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
	trait_df <- read.table(trait_gene_pip_summary_file, header=TRUE,sep="\t")

	#ordered_tissues <- sort(unique(as.character(trait_df$tissue_name)))

	tmper <- c()
	for (pip_iter in 1:length(pip_threshs)) {
		pip_thresh <- pip_threshs[pip_iter]
		tmp_df <- trait_df[trait_df$PIP > pip_thresh,]
		#print(head(tmp_df))
		for (tissue_iter in 1:length(ordered_tissues)) {
			tissue_name <- ordered_tissues[tissue_iter]
			pip_sum = length(tmp_df$PIP[as.character(tmp_df$tissue_name) ==tissue_name])
			tissue_vec <- c(tissue_vec, tissue_name)
			count_vec <- c(count_vec, pip_sum)
			pip_thresh_vec <- c(pip_thresh_vec, pip_thresh)
			if (pip_thresh == pip_threshs[1]) {
				tmper <- c(tmper, pip_sum)
			}
		}


	}

	ord <- order(-tmper)
	df <- data.frame(expected_causal_genes=count_vec, tissue=factor(tissue_vec,levels=ordered_tissues[ord]), pip_thresh=factor(pip_thresh_vec))
	df$tissue = str_replace_all(as.character(df$tissue), "-", "_")
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="adipose sub.", Adipose_Visceral_Omentum="adipose visceral", Breast_Mammary_Tissue="breast", Cells_Cultured_fibroblasts="fibroblast",Heart_Atrial_Appendage="heart atrial",Skin_Sun_Exposed_Lower_leg="skin (sun)",Skin_Not_Sun_Exposed_Suprapubic="skin (no sun)", Small_Intestine_Terminal_Ileum="small intestine", Brain_Anterior_cingulate_cortex_BA24="brain anterior cortex", Brain_Nucleus_accumbens_basal_ganglia="brain basal ganglia", Esophagus_Gastroesophageal_Junction="esophagus gastro. junct.", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="brain spinal chord", Whole_Blood="whole blood", Adrenal_Gland="adrenal gland", Spleen="spleen", Brain_BasalGanglia="brain basal ganglia", Thyroid="thyroid", Artery_Aorta="artery aorta", Artery_Coronary="artery coronary", Artery_Tibial="artery tibial", Brain_Cerebellum="brain cerebellum", Brain_Cortex="brain cortex", Brain_Limbic="brain limbic", Brain_Substantia_nigra="brain subst. nigra", Colon_Sigmoid="colon sigmoid", Colon_Transverse="colon transverse", Esophagus_Mucosa="esophagus mucosa", Esophagus_Muscularis="esophagus muscularis", Heart_Left_Ventricle="heart left ventricle", Liver="liver", Lung="lung", Minor_Salivary_Gland="minor salivary gland", Muscle_Skeletal="muscle skeletal", Nerve_Tibial="nerve tibial", Ovary="ovary", Prostate="prostate", Stomach="stomach", Uterus="uterus", Vagina="vagina")
	ordered_tissues2 <- as.character(df$tissue)[1:length(unique(df$tissue))]
	df$tissue = factor(df$tissue, levels=ordered_tissues2[ord])

	red_color=brewer.pal(n = 9, name = "Reds")[7]

	p<-ggplot(df, aes(x=tissue, y=expected_causal_genes)) +
  		geom_bar(stat="identity",position=position_dodge(), fill=red_color)+figure_theme() +
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  		theme(legend.position="bottom") +
  		labs(x="", y="No. gene-tissue pairs\n(PIP > 0.5)",title=trait_name_readable)

  	return(p)
}


make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait <- function(trait_names, trait_names_readable, method_version, tgfm_results_dir, gene_type) {
	# First extract data from each trait
	trait_vec <- c()
	tissue_vec <- c()
	broad_tissue_vec <- c()
	pip_vec <- c()
	# Extract data from each trait
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]

		trait_gene_pip_summary_file <- paste0(tgfm_results_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_tgfm_per_gene_tissue_pip_summary.txt")
		print(trait_gene_pip_summary_file)
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


make_number_of_high_pip_gene_tissue_pairs_method_comparison_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits) {

	trait_names_vec <- c()
	number_elements_vec <- c()
	element_type_vec <- c()



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {


		gene_type <- "component_gene"
		summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_threshold, ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		gene_tissue_count = tmp$count[2]
		
		gene_type <- "all_non_zero_gene"
		summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_threshold, ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		gene_count = tmp$count[2]

		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		number_elements_vec <- c(number_elements_vec, gene_tissue_count)
		number_elements_vec <- c(number_elements_vec, gene_count)
		element_type_vec <- c(element_type_vec, "component_gene")
		element_type_vec <- c(element_type_vec, "all_non_zero_gene")
		}
	}

	df <- data.frame(trait=trait_names_vec, number_of_elements=number_elements_vec, element_type=factor(element_type_vec, levels=c("component_gene", "all_non_zero_gene")))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$element_type) == "component_gene",]

	indices = order(-tmp_df$number_of_elements)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	#red_color=brewer.pal(n = 9, name = "Reds")[7]


	p <- ggplot(df, aes(x=trait, y=number_of_elements, fill=element_type)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="No. fine-mapped\ngenetic elements", x="", fill="", title=paste0("PIP >= ", pip_threshold)) +
    		figure_theme() +
    		theme(legend.position=c(0.8, 0.86))
    return(p)


}

make_number_of_high_pip_genes_vs_gene_tissues_barplot <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits) {

	trait_names_vec <- c()
	number_elements_vec <- c()
	element_type_vec <- c()



	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {


		summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genetic_elements_pip_", pip_threshold, ".txt")
		tmp <- read.table(summary_file, header=TRUE)
		gene_tissue_count = tmp$count[2]
		
		summary_file2 <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genes_pip_", pip_threshold, ".txt")
		tmp <- read.table(summary_file2, header=TRUE)
		gene_count = tmp$count[1]

		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		trait_names_vec <- c(trait_names_vec, trait_name_readable)
		number_elements_vec <- c(number_elements_vec, gene_tissue_count)
		number_elements_vec <- c(number_elements_vec, gene_count)
		element_type_vec <- c(element_type_vec, "Gene-Tissue")
		element_type_vec <- c(element_type_vec, "Gene")
		}
	}

	df <- data.frame(trait=trait_names_vec, number_of_elements=number_elements_vec, element_type=factor(element_type_vec, levels=c("Gene", "Gene-Tissue")))
	df$trait <- recode(df$trait, biochemistry_Cholesterol="Cholesterol", biochemistry_LDLdirect="LDLdirect", biochemistry_VitaminD="VitaminD", blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT="Reticulocyte count", blood_MEAN_CORPUSCULAR_HEMOGLOBIN="Hemoglobin", blood_MEAN_PLATELET_VOL="Platelet volume", blood_MONOCYTE_COUNT="Monocyte count", body_BMIz="BMI", body_HEIGHTz="Height", body_WHRadjBMIz="WHR-adj BMI", bp_DIASTOLICadjMEDz="Diastolic BP", cov_EDU_COLLEGE="College Education", disease_ALLERGY_ECZEMA_DIAGNOSED="Eczema", disease_CARDIOVASCULAR="Cardiovascular", lung_FEV1FVCzSMOKE="FEV1FVCz", repro_MENARCHE_AGE="Menarche Age")

	tmp_df = df[as.character(df$element_type) == "Gene",]

	indices = order(-tmp_df$number_of_elements)
	df$trait = factor(df$trait, levels=as.character(tmp_df$trait)[indices])

	red_color=brewer.pal(n = 9, name = "Reds")[7]

	p <- ggplot(df, aes(x=trait, y=number_of_elements, fill=element_type)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		#geom_bar(stat="identity", position="stack") +
    		scale_fill_manual(values=c("palegreen4", red_color)) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=11)) +
    		labs(y="No. fine-mapped\ngenetic elements", x="", fill="", title=paste0("PIP >= ", pip_threshold)) +
    		figure_theme() +
    		theme(legend.position=c(0.8, 0.86))
    return(p)


}

make_pip_density_histogram_across_element_classes <- function(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits) {
	pip_vec <- c()
	element_class_vec <- c()
	trait_vec <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_name_readable <- trait_names_readable[trait_iter]
		if (trait_name %in% independent_traits) {

			# Add variants
			variant_summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_variants_cross_pip_threshold_sqrt_plot_input.txt")
			tmp_data <- read.table(variant_summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="variant",]
			tmp_data = tmp_data[tmp_data$PIP_threshold >= .2,]
			pip_vec <- c(pip_vec, tmp_data$PIP_threshold)
			n_variants = length(tmp_data$PIP_threshold)
			element_class_vec <- c(element_class_vec, rep("Variant", n_variants))
			trait_vec <- c(trait_vec, rep(trait_name_readable, n_variants))

			# Add genes
			gene_summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_genes_cross_pip_threshold_sqrt_plot_input.txt")
			tmp_data <- read.table(gene_summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="gene",]
			tmp_data = tmp_data[tmp_data$PIP_threshold >= .2,]
			pip_vec <- c(pip_vec, tmp_data$PIP_threshold)
			n_genes = length(tmp_data$PIP_threshold)
			element_class_vec <- c(element_class_vec, rep("Gene", n_genes))
			trait_vec <- c(trait_vec, rep(trait_name_readable, n_genes))

			# Add gene-tissue pairs
			gene_tissue_summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_gene_tissue_pairs_cross_pip_threshold_sqrt_plot_input.txt")
			tmp_data <- read.table(gene_tissue_summary_file, header=TRUE, sep="\t")
			tmp_data <- tmp_data[as.character(tmp_data$element_class)=="gene",]
			tmp_data = tmp_data[tmp_data$PIP_threshold >= .2,]
			pip_vec <- c(pip_vec, tmp_data$PIP_threshold)
			n_genes = length(tmp_data$PIP_threshold)
			element_class_vec <- c(element_class_vec, rep("Gene-Tissue", n_genes))
			trait_vec <- c(trait_vec, rep(trait_name_readable, n_genes))
		}
	}

	df <- data.frame(PIP=pip_vec, genetic_element=factor(element_class_vec, levels=c("Gene-Tissue", "Gene", "Variant")), trait=trait_vec)

	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[7]
	red_color=brewer.pal(n = 9, name = "Reds")[7]


	pp <- ggplot(df, aes(x = PIP, color=genetic_element)) +
		geom_density(size=1.5) +
		figure_theme() +
		labs(x="TGFM PIP", color="") +
		scale_colour_manual(values=c(red_color, green_color, blue_color))
	return(pp)
}


mean_se_barplot_of_pops_score_comparison_binned_by_tgfm_pip <- function(df_full_1, df_full_2, method_1_name, method_2_name, independent_traits) {

	indices = df_full_1$trait_name %in% independent_traits
	df_1 = df_full_1[indices,]
	indices = df_full_2$trait_name %in% independent_traits
	df_2 = df_full_2[indices,]

	pops_mean_vec <- c()
	pops_mean_se_vec <- c()
	bin_names_vec <- c()
	method_vec <- c()

	# METHOD 1
	threshold_lb <- 0.0
	threshold_ub <- .01
	indices = (df_1$tgfm_gene_pip >= threshold_lb) & (df_1$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_1$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_1$pops_score[indices])/sqrt(length(df_1$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_1_name)

	threshold_lb <- .01
	threshold_ub <- .25
	indices = (df_1$tgfm_gene_pip >= threshold_lb) & (df_1$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_1$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_1$pops_score[indices])/sqrt(length(df_1$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_1_name)


	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df_1$tgfm_gene_pip >= threshold_lb) & (df_1$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_1$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_1$pops_score[indices])/sqrt(length(df_1$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_1_name)

	threshold_lb <- .5
	threshold_ub <- .7
	indices = (df_1$tgfm_gene_pip >= threshold_lb) & (df_1$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_1$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_1$pops_score[indices])/sqrt(length(df_1$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_1_name)

	threshold_lb <- .7
	threshold_ub <- .9
	indices = (df_1$tgfm_gene_pip >= threshold_lb) & (df_1$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_1$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_1$pops_score[indices])/sqrt(length(df_1$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_1_name)


	threshold_lb <- .9
	threshold_ub <- 1.0
	indices = (df_1$tgfm_gene_pip >= threshold_lb) & (df_1$tgfm_gene_pip <= threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_1$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_1$pops_score[indices])/sqrt(length(df_1$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP <= ", threshold_ub))
	method_vec <- c(method_vec, method_1_name)


	# METHOD 2
	threshold_lb <- 0.0
	threshold_ub <- .01
	indices = (df_2$tgfm_gene_pip >= threshold_lb) & (df_2$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_2$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_2$pops_score[indices])/sqrt(length(df_2$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_2_name)

	threshold_lb <- .01
	threshold_ub <- .25
	indices = (df_2$tgfm_gene_pip >= threshold_lb) & (df_2$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_2$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_2$pops_score[indices])/sqrt(length(df_2$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_2_name)


	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df_2$tgfm_gene_pip >= threshold_lb) & (df_2$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_2$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_2$pops_score[indices])/sqrt(length(df_2$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_2_name)

	threshold_lb <- .5
	threshold_ub <- .7
	indices = (df_2$tgfm_gene_pip >= threshold_lb) & (df_2$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_2$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_2$pops_score[indices])/sqrt(length(df_2$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_2_name)

	threshold_lb <- .7
	threshold_ub <- .9
	indices = (df_2$tgfm_gene_pip >= threshold_lb) & (df_2$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_2$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_2$pops_score[indices])/sqrt(length(df_2$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))
	method_vec <- c(method_vec, method_2_name)


	threshold_lb <- .9
	threshold_ub <- 1.0
	indices = (df_2$tgfm_gene_pip >= threshold_lb) & (df_2$tgfm_gene_pip <= threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df_2$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df_2$pops_score[indices])/sqrt(length(df_2$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP <= ", threshold_ub))
	method_vec <- c(method_vec, method_2_name)



	df2 <- data.frame(pops=pops_mean_vec, pops_se=pops_mean_se_vec, tgfm_bin=factor(bin_names_vec), method=factor(method_vec, levels=c("component_gene", "all_non_zero_gene")))
	print(df2)

	green_color=brewer.pal(n = 9, name = "Greens")[6]

	p <- ggplot(df2) +
    		geom_bar( aes(x=tgfm_bin, y=pops, fill=method), stat="identity", position="dodge") +
    		labs(y="PoPS score", x="") +
    		geom_errorbar( aes(x=tgfm_bin,fill=method, ymin=pops-(1.96*pops_se), ymax=pops+(1.96*pops_se)), width=0.4, colour="grey45", alpha=0.9, size=1.0,position=position_dodge(.9)) +
    		theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    		theme(legend.position="bottom")+
    		figure_theme()



    return(p)


}


mean_se_barplot_of_pops_score_binned_by_tgfm_pip <- function(df_full, independent_traits) {

	indices = df_full$trait_name %in% independent_traits
	df = df_full[indices,]


	pops_mean_vec <- c()
	pops_mean_se_vec <- c()
	bin_names_vec <- c()

	threshold_lb <- 0.0
	threshold_ub <- .01
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))

	threshold_lb <- .01
	threshold_ub <- .25
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))


	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))

	threshold_lb <- .5
	threshold_ub <- .7
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))

	threshold_lb <- .7
	threshold_ub <- .9
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))


	threshold_lb <- .9
	threshold_ub <- 1.0
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip <= threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP <= ", threshold_ub))


	df2 <- data.frame(pops=pops_mean_vec, pops_se=pops_mean_se_vec, tgfm_bin=factor(bin_names_vec))

	green_color=brewer.pal(n = 9, name = "Greens")[6]

	p <- ggplot(df2) +
    		geom_bar( aes(x=tgfm_bin, y=pops), stat="identity", fill=green_color, alpha=.95) +
    		#theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="PoPS score", x="") +
    		geom_errorbar( aes(x=tgfm_bin, ymin=pops-(1.96*pops_se), ymax=pops+(1.96*pops_se)), width=0.4, colour="grey45", alpha=0.9, size=1.0) +
    		theme(axis.text.x = element_text(angle = 45, hjust=1,size=11)) + 
    		figure_theme()

    return(p)
}

mean_se_barplot_of_pops_score_binned_by_tgfm_pip_supp_table <- function(df_full, independent_traits) {

	indices = df_full$trait_name %in% independent_traits
	df = df_full[indices,]


	pops_mean_vec <- c()
	pops_mean_se_vec <- c()
	bin_names_vec <- c()

	threshold_lb <- 0.0
	threshold_ub <- .01
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))

	threshold_lb <- .01
	threshold_ub <- .25
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))


	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))

	threshold_lb <- .5
	threshold_ub <- .7
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))

	threshold_lb <- .7
	threshold_ub <- .9
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP < ", threshold_ub))


	threshold_lb <- .9
	threshold_ub <- 1.0
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip <= threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= PIP <= ", threshold_ub))


	df2 <- data.frame(tgfm_bin=factor(bin_names_vec), pops=pops_mean_vec, pops_se=pops_mean_se_vec)




    return(df2)


}


count_up_number_of_genetic_elements <- function(trait_names,tgfm_organized_results_dir, pip_threshold) {
	gts = 0
	genes = 0
	variants = 0
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		summary_file1 <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_n_causal_genetic_elements_pip_", pip_threshold, ".txt")
		summary_file2 <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_n_causal_genes_pip_", pip_threshold, ".txt")

		aa <- read.table(summary_file1, header=TRUE)
		variants = variants + aa$count[1]
		gts = gts + aa$count[2]
		bb = read.table(summary_file2, header=TRUE)
		genes = genes + bb$count[1]
	}

	print(paste0("Gene-Tissue count: ", gts))
	print(paste0("Gene count: ", genes))
	print(paste0("Variant count: ", variants))
}



make_non_disease_specific_gene_set_enrichment_barplot_cross_pip_thresholds <- function(non_disease_specific_gene_set_enrichment_dir) {
	# Enrichment summary file
	enrichment_summary_file = paste0(non_disease_specific_gene_set_enrichment_dir, "global_enrichment_summary.txt")
	df <- read.table(enrichment_summary_file, header=TRUE, sep="\t")
	
	df2 = df[df$PIP==.5,]


    ordering = order(-df2$odds_ratio)
    df$gene_set_name = factor(df$gene_set_name, levels=as.character(df2$gene_set_name)[ordering])

    df$PIP = factor(df$PIP)

	red_color=brewer.pal(n = 9, name = "Greens")[7]
	red_color1=brewer.pal(n = 9, name = "Greens")[3]
	red_color2=brewer.pal(n = 9, name = "Greens")[5]


	p <- ggplot(df) +
    		geom_bar( aes(x=gene_set_name, y=odds_ratio, fill=PIP), stat="identity", position="dodge") +
    		#theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="Odds ratio", x="") +
    		geom_errorbar( aes(x=gene_set_name,fill=PIP, ymin=odds_ratio_lb, ymax=odds_ratio_ub), width=0.4, colour="grey45", alpha=0.9, size=1.0,position=position_dodge(.9)) +
    		theme(axis.text.x = element_text(angle = 64, hjust=1)) + 
    		scale_fill_manual(values=c(red_color1, red_color2, red_color))+
    		geom_hline(yintercept=1) + 
    		figure_theme()


    return(p)
}

make_non_disease_specific_gene_set_enrichment_suppData_cross_pip_thresh <- function(non_disease_specific_gene_set_enrichment_dir) {
	# Enrichment summary file
	enrichment_summary_file = paste0(non_disease_specific_gene_set_enrichment_dir, "global_enrichment_summary.txt")
	df <- read.table(enrichment_summary_file, header=TRUE, sep="\t")
	#df = df[df$PIP==.5,]

	#green_color=brewer.pal(n = 9, name = "Greens")[6]


    #ordering = order(-df$odds_ratio)
    #df$gene_set_name = factor(df$gene_set_name, levels=as.character(df$gene_set_name)[ordering])

    df2 = data.frame(gene_set_name=df$gene_set_name, PIP=df$PIP,odds_ratio=df$odds_ratio, odds_ratio_95_ci_lb=df$odds_ratio_lb, odds_ratio_95_ci_ub=df$odds_ratio_ub, odds_ratio_pvalue=df$odds_ratio_pvalue)

    #df2$gene_set_name = factor(df2$gene_set_name, levels=as.character(df$gene_set_name)[ordering])

    return(df2)
}


make_non_disease_specific_gene_set_enrichment_barplot_at_single_pip_thresh <- function(non_disease_specific_gene_set_enrichment_dir) {
	# Enrichment summary file
	enrichment_summary_file = paste0(non_disease_specific_gene_set_enrichment_dir, "global_enrichment_summary.txt")
	df <- read.table(enrichment_summary_file, header=TRUE, sep="\t")
	df = df[df$PIP==.5,]

	green_color=brewer.pal(n = 9, name = "Greens")[6]


    ordering = order(-df$odds_ratio)
    df$gene_set_name = factor(df$gene_set_name, levels=as.character(df$gene_set_name)[ordering])

	p <- ggplot(df) +
    		geom_bar( aes(x=gene_set_name, y=odds_ratio), stat="identity", fill=green_color, alpha=.95) +
    		#theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="Odds ratio", x="") +
    		geom_errorbar( aes(x=gene_set_name, ymin=odds_ratio_lb, ymax=odds_ratio_ub), width=0.4, colour="grey45", alpha=0.9, size=1.0) +
    		theme(axis.text.x = element_text(angle = 64, hjust=1)) + 
    		geom_hline(yintercept=1) + 
    		figure_theme()


    return(p)
}

make_non_disease_specific_gene_set_enrichment_method_comparison_barplot_at_single_pip_thresh <- function(non_disease_specific_gene_set_enrichment_dir, pip_thresh) {

	# Enrichment summary file
	enrichment_summary_file = paste0(non_disease_specific_gene_set_enrichment_dir, "global_enrichment_summary.txt")
	df1 <- read.table(enrichment_summary_file, header=TRUE, sep="\t")
	df1 = df1[df1$PIP==pip_thresh,]

	enrichment_summary_file = paste0(non_disease_specific_gene_set_enrichment_dir, "all_non_zero_gene_global_enrichment_summary.txt")
	df2 <- read.table(enrichment_summary_file, header=TRUE, sep="\t")
	df2 = df2[df2$PIP==pip_thresh,]



	gene_set_name = c(as.character(df1$gene_set_name), as.character(df2$gene_set_name))
	odds_ratio <- c(df1$odds_ratio, df2$odds_ratio)
	odds_ratio_lb <- c(df1$odds_ratio_lb, df2$odds_ratio_lb)
	odds_ratio_ub <- c(df1$odds_ratio_ub, df2$odds_ratio_ub)
	method <- c(rep("component_gene", length(df1$gene_set_name)), rep("all_non_zero_gene", length(df2$gene_set_name)))

	df <- data.frame(gene_set_name=gene_set_name, odds_ratio=odds_ratio, odds_ratio_ub=odds_ratio_ub, odds_ratio_lb=odds_ratio_lb, method=factor(method,levels=c("component_gene", "all_non_zero_gene")))


    ordering = order(-df1$odds_ratio)
    df$gene_set_name = factor(as.character(df$gene_set_name), levels=as.character(df$gene_set_name)[ordering])

	p2 <- ggplot(df) +
    		geom_bar( aes(x=gene_set_name, y=odds_ratio, fill=method), stat="identity", position="dodge") +
    		#theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="Odds ratio", x="",title=paste0("PIP=", pip_thresh)) +
    		geom_errorbar( aes(x=gene_set_name, fill=method, ymin=odds_ratio_lb, ymax=odds_ratio_ub), width=0.4, colour="grey45", alpha=0.9, size=1.0,position=position_dodge(.9)) +
    		theme(axis.text.x = element_text(angle = 64, hjust=1)) + 
    		geom_hline(yintercept=1) + 
    		figure_theme()






    return(p2)
}





make_tissue_tissue_correlation_sample_size_scatterplot <- function(tissue_tissue_correlation_file, tissue_names, tissue_ss) {
	df <- read.table(tissue_tissue_correlation_file, header=TRUE)
	nrow = dim(df)[1]

	tissue_names <- as.character(tissue_names)
	indices = tissue_names != "Testis"
	tissue_names = tissue_names[indices]
	tissue_ss = tissue_ss[indices]
	avg_corr_vec <- c()

	for (tiss_iter in 1:length(tissue_names)) {
		tissue_name = tissue_names[tiss_iter]

		tmp_df = df[df$tissue_ii == tissue_name,]
		tmp_df = tmp_df[tmp_df$tissue_jj != tissue_name,]

		mean_corry = mean(tmp_df$average_correlation)
		avg_corr_vec <- c(avg_corr_vec, mean_corry)
	}

	df2 <- data.frame(tissue=tissue_names, sample_size=tissue_ss, average_correlation=avg_corr_vec)



	pp <- ggplot(df2, aes(x=sample_size, y=average_correlation)) + geom_point() +
		figure_theme() +
		labs(x="Tissue sample size", y="Average correlation")
	return(pp)
}


make_number_of_high_pip_sc_gene_tissue_pairs_in_each_trait_heatmap_barplot <- function(trait_name, trait_name_readable, method_version, tgfm_organized_results_dir, cell_type_names, upper_bound=25) {
	cell_type_names = as.character(cell_type_names)
	tissue_names_vec <- c()
	number_elements_vec <- c()
	pip_threshold_vec <- c()

	total_hits_vec <- c()
	mid_point_vec <- c()
	tissue_names_vec2 <- c()

	for (cell_type_iter in 1:length(cell_type_names)) {
		cell_type_name <- cell_type_names[cell_type_iter]


		summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_", method_version, "_tgfm_n_causal_tissue_gene_tissue_pairs_", cell_type_name, "_cross_pip_threshold_sqrt_plot_input.txt")
		tmp_data <- read.table(summary_file, header=TRUE, sep="\t")
		tmp_data <- tmp_data[as.character(tmp_data$element_class)=="gene",]

		tmp_data = tmp_data[tmp_data$PIP_threshold >= .2,]
			
		number_elements_vec <- c(number_elements_vec, tmp_data$n_elements)
		pip_threshold_vec <- c(pip_threshold_vec, tmp_data$PIP_threshold)
		tissue_names_vec <- c(tissue_names_vec, rep(cell_type_name, length(tmp_data$PIP_threshold)))

		total_hits_vec <- c(total_hits_vec, sum(tmp_data$n_elements))
		tissue_names_vec2 <- c(tissue_names_vec2, cell_type_name)
		mid_point_vec <- c(mid_point_vec, sum(tmp_data$n_elements[tmp_data$PIP_threshold >= .5]))


	}

	df <- data.frame(trait=tissue_names_vec, PIP=pip_threshold_vec, num_elements=number_elements_vec)

	indices = order(-mid_point_vec, -total_hits_vec)


	df$trait = str_replace_all(as.character(df$trait), "-", "_")
	df$trait <- recode(df$trait, Adipose_Subcutaneous="adipose sub.",Pituitary="pituitary", Adipose_Visceral_Omentum="adipose visceral", Breast_Mammary_Tissue="breast", Cells_Cultured_fibroblasts="fibroblast",Heart_Atrial_Appendage="heart atrial",Skin_Sun_Exposed_Lower_leg="skin (sun)",Skin_Not_Sun_Exposed_Suprapubic="skin (no sun)", Small_Intestine_Terminal_Ileum="small intestine", Brain_Anterior_cingulate_cortex_BA24="brain anterior cortex", Brain_Nucleus_accumbens_basal_ganglia="brain basal ganglia", Esophagus_Gastroesophageal_Junction="esophagus gastro. junct.", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="brain spinal chord", Whole_Blood="whole blood", Adrenal_Gland="adrenal gland", Spleen="spleen", Brain_BasalGanglia="brain basal ganglia", Thyroid="thyroid", Artery_Aorta="artery aorta", Artery_Coronary="artery coronary", Artery_Tibial="artery tibial", Brain_Cerebellum="brain cerebellum", Brain_Cortex="brain cortex", Brain_Limbic="brain limbic", Brain_Substantia_nigra="brain subst. nigra", Colon_Sigmoid="colon sigmoid", Colon_Transverse="colon transverse", Esophagus_Mucosa="esophagus mucosa", Esophagus_Muscularis="esophagus muscularis", Heart_Left_Ventricle="heart left ventricle", Liver="liver", Lung="lung", Minor_Salivary_Gland="minor salivary gland", Muscle_Skeletal="muscle skeletal", Nerve_Tibial="nerve tibial", Ovary="ovary", Prostate="prostate", Stomach="stomach", Uterus="uterus", Vagina="vagina")

	tissue_names_vec2 = str_replace_all(as.character(tissue_names_vec2), "-", "_")
	tissue_names_vec2 = recode(tissue_names_vec2, Adipose_Subcutaneous="adipose sub.",Pituitary="pituitary", Adipose_Visceral_Omentum="adipose visceral", Breast_Mammary_Tissue="breast", Cells_Cultured_fibroblasts="fibroblast",Heart_Atrial_Appendage="heart atrial",Skin_Sun_Exposed_Lower_leg="skin (sun)",Skin_Not_Sun_Exposed_Suprapubic="skin (no sun)", Small_Intestine_Terminal_Ileum="small intestine", Brain_Anterior_cingulate_cortex_BA24="brain anterior cortex", Brain_Nucleus_accumbens_basal_ganglia="brain basal ganglia", Esophagus_Gastroesophageal_Junction="esophagus gastro. junct.", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="brain spinal chord", Whole_Blood="whole blood", Adrenal_Gland="adrenal gland", Spleen="spleen", Brain_BasalGanglia="brain basal ganglia", Thyroid="thyroid", Artery_Aorta="artery aorta", Artery_Coronary="artery coronary", Artery_Tibial="artery tibial", Brain_Cerebellum="brain cerebellum", Brain_Cortex="brain cortex", Brain_Limbic="brain limbic", Brain_Substantia_nigra="brain subst. nigra", Colon_Sigmoid="colon sigmoid", Colon_Transverse="colon transverse", Esophagus_Mucosa="esophagus mucosa", Esophagus_Muscularis="esophagus muscularis", Heart_Left_Ventricle="heart left ventricle", Liver="liver", Lung="lung", Minor_Salivary_Gland="minor salivary gland", Muscle_Skeletal="muscle skeletal", Nerve_Tibial="nerve tibial", Ovary="ovary", Prostate="prostate", Stomach="stomach", Uterus="uterus", Vagina="vagina")

	df$trait = factor(df$trait, levels=as.character(tissue_names_vec2)[indices])
	df2 <- data.frame(trait=factor(as.character(tissue_names_vec2), levels=as.character(tissue_names_vec2)[indices]), midpoint=mid_point_vec)
	df2$trait = factor(df2$trait, levels=as.character(tissue_names_vec2)[indices])

	valid_tissues = as.character(tissue_names_vec2)[indices][1:10]

	df = df[df$trait %in% valid_tissues,]

	df2 = df2[df2$trait %in% valid_tissues,]

	#ordered_tissues2 <- as.character(df$trait)[1:length(unique(df$trait))]
	#df$trait = factor(df$trait, levels=trait[ord])




   p <- ggplot() +
  	geom_bar(data=df, aes(x = trait,y = num_elements,fill = PIP, color=PIP), stat="identity", width=.9) + 
  	figure_theme() +
  	#scale_y_continuous(breaks=c(0.0,sqrt(5), sqrt(20), sqrt(50), sqrt(100), sqrt(200), sqrt(400), sqrt(600)), labels=c(0, 5, 20, 50, 100, 200, 400,600)) +
  	scale_y_continuous(breaks=c(0.0, sqrt(10), sqrt(50), sqrt(100)), labels=c("0", "10", "50", "100")) +
  	theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size=11)) +
  	scale_fill_distiller(direction=1, palette = "Purples", limits=c(.2,1.0)) +
  	scale_color_distiller(direction=1, palette = "Purples", limits=c(.2,1.0)) +	
  	theme(plot.title = element_text(hjust = 0.5)) +
  	geom_errorbar(data=df2, aes(x=trait,y=midpoint, ymax=midpoint, ymin=midpoint)) +
  	 labs(x="", y="No. fine-mapped\ngene-tissue pairs", title=trait_name_readable)
 	return(p)


}
make_trait_tissue_chromatin_overlap_density <- function(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE) {
	# Load in data
	df <- read.table(trait_tissue_chromatin_overlap_barplot, header=TRUE, sep="\t")


	df$sldsc_p = -log10(df$sldsc_p)

	# Potentially filter to independent traits
	if (independent_traits) {
		df <- df[as.character(df$trait_type) == "independent", ]
	}

	arr = rep("FDR > 0.05", length(df$tgfm_sig))
	arr[as.character(df$tgfm_sig) == "True"] = "FDR <= 0.05"
	df$tgfm_fdr = arr
	df$tgfm_fdr = factor(df$tgfm_fdr, levels=c("FDR <= 0.05", "FDR > 0.05"))
	red_color =brewer.pal(n = 9, name = "Reds")[6]

	pp <- ggplot(df, aes(x=sldsc_p,  fill=tgfm_fdr,colour = tgfm_fdr)) +
  		geom_density(alpha = 0.2) +
    	scale_fill_manual(values=c(red_color, "black")) +
  		scale_colour_manual(values=c(red_color, "grey")) +
  		figure_theme() + 
    	theme(legend.position="top") + 
  		labs(x="SLDSC -log10(p-value)",fill="Tissue-trait pairs (TGFM tissue-specific prior)",colour="Tissue-trait pairs (TGFM tissue-specific prior)")

  	return(pp)

}
get_trait_tissue_chromatin_overlap_cdf_raw_data <- function(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE) {
	# Load in data
	df <- read.table(trait_tissue_chromatin_overlap_barplot, header=TRUE, sep="\t")


	df$sldsc_p = -log10(df$sldsc_p)

	# Potentially filter to independent traits
	if (independent_traits) {
		df <- df[as.character(df$trait_type) == "independent", ]
	}

	arr = rep("FDR > 0.05", length(df$tgfm_sig))
	arr[as.character(df$tgfm_sig) == "True"] = "FDR <= 0.05"
	df$tgfm_fdr = arr
	df$tgfm_fdr = factor(df$tgfm_fdr, levels=c("FDR <= 0.05", "FDR > 0.05"))

  	return(df)

}

make_trait_tissue_chromatin_overlap_cdf <- function(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE) {
	# Load in data
	df <- read.table(trait_tissue_chromatin_overlap_barplot, header=TRUE, sep="\t")


	df$sldsc_p = -log10(df$sldsc_p)

	# Potentially filter to independent traits
	if (independent_traits) {
		df <- df[as.character(df$trait_type) == "independent", ]
	}

	null_stringer <- paste0("FDR > ", fdr)
	alt_stringer <- paste0("FDR <= ", fdr)


	arr = rep(null_stringer, length(df$tgfm_sig))
	arr[as.character(df$tgfm_sig) == "True"] = alt_stringer
	df$tgfm_fdr = arr
	df$tgfm_fdr = factor(df$tgfm_fdr, levels=c(alt_stringer, null_stringer))
	red_color =brewer.pal(n = 9, name = "Reds")[6]

	pp <- ggplot(df, aes(x=sldsc_p,  fill=tgfm_fdr,colour = tgfm_fdr, y= 1-..y..)) +
  		stat_ecdf(size=1.4) +
  		figure_theme() + 
     	scale_fill_manual(values=c(red_color, "black")) +
  		scale_colour_manual(values=c(red_color, "grey37")) +
  		theme(legend.position="top") + 
  		geom_vline(xintercept = -log10(.05),linetype="dotted") +
  		labs(x="SLDSC -log10(p-value) threshold",y="Proportion SLDSC significant", fill="Tissue-trait pairs (TGFM tissue-specific prior)",colour="Tissue-trait pairs (TGFM tissue-specific prior)" )

  	return(pp)

}


make_trait_tissue_chromatin_overlap_histogram <- function(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE) {
	# Load in data
	df <- read.table(trait_tissue_chromatin_overlap_barplot, header=TRUE, sep="\t")

	df$sldsc_p = -log10(df$sldsc_p)

	# Potentially filter to independent traits
	if (independent_traits) {
		df <- df[as.character(df$trait_type) == "independent", ]
	}

	
	arr = rep("FDR > 0.05", length(df$tgfm_sig))
	arr[as.character(df$tgfm_sig) == "True"] = "FDR <= 0.05"
	df$tgfm_fdr = arr
	df$tgfm_fdr = factor(df$tgfm_fdr, levels=c("FDR <= 0.05", "FDR > 0.05"))
	red_color =brewer.pal(n = 9, name = "Reds")[6]

	pp <- ggplot(df, aes(x=sldsc_p, colour = tgfm_fdr,fill=tgfm_fdr)) +
  		geom_histogram(aes(y=..density..),alpha = 0.2, position='identity') +
  		figure_theme() + 
  		scale_fill_manual(values=c(red_color, "black")) +
  		scale_colour_manual(values=c(red_color, "grey")) +
  		theme(legend.position="top") +
  		labs(x="SLDSC -log10(p-value)", colour="Tissue-trait pairs (TGFM tissue-specific prior)", fill="Tissue-trait pairs (TGFM tissue-specific prior)") 
  	return(pp)


}

make_trait_tissue_chromatin_overlap_barplot <- function(trait_tissue_chromatin_overlap_barplot, independent_traits=TRUE) {
	# Load in data
	df <- read.table(trait_tissue_chromatin_overlap_barplot, header=TRUE, sep="\t")

	# Potentially filter to independent traits
	if (independent_traits) {
		df <- df[as.character(df$trait_type) == "independent", ]
	}

	# Now split up data frame into trait-tissue pairs identified as significant by TGFM
	tgfm_sig_df <- df[as.character(df$tgfm_sig) == "True",]
	tgfm_null_df <- df[as.character(df$tgfm_sig) == "False",]

	aa = sum((tgfm_sig_df$sldsc_p < 0.05) & (tgfm_sig_df$sldsc_z > 0.0))
	bb = dim(tgfm_sig_df)[1] - aa
	n_tgfm_sig = aa +bb

	cc = sum((tgfm_null_df$sldsc_p < 0.05) & (tgfm_null_df$sldsc_z > 0.0))
	dd = dim(tgfm_null_df)[1] - cc
	n_tgfm_null = cc + dd

	plot_df <- data.frame(x=c("FDR <= 0.05", "FDR > 0.05"), y=c(aa/n_tgfm_sig, cc/n_tgfm_null), totals=c(paste0(n_tgfm_sig), paste0(n_tgfm_null)))

	red_color =brewer.pal(n = 9, name = "Reds")[6]


	pp <- ggplot(data=plot_df, aes(x=x, y=y, fill=x)) +
  		geom_bar(stat="identity", width=0.9) +
  		geom_text(aes(label=totals, y=c(0.033,0.033)),vjust=0.0) +
  		labs(x="Tissue-trait pairs\n(TGFM tissue-specific prior)", y="Replication\n rate in SLDSC") +
  		figure_theme() +
  		scale_fill_manual(values=c(red_color,red_color)) + 
  		theme(legend.position="None") +
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) 
  	return(pp)
}

make_trait_tissue_chromatin_overlap_barplot_multiple_thresh <- function(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=TRUE) {
	# Load in data
	df <- read.table(trait_tissue_chromatin_overlap_barplot, header=TRUE, sep="\t")

	# Potentially filter to independent traits
	if (independent_traits) {
		df <- df[as.character(df$trait_type) == "independent", ]
	}


	fdr_vec <- c()
	replication_frac_vec <- c()
	totals_vec <- c()
	pval_vec <- c()

	# Now split up data frame into trait-tissue pairs identified as significant by TGFM
	tgfm_sig_df <- df[as.character(df$tgfm_sig) == "True",]
	tgfm_null_df <- df[as.character(df$tgfm_sig) == "False",]

	aa = sum((tgfm_sig_df$sldsc_p < 0.05) & (tgfm_sig_df$sldsc_z > 0.0))
	bb = dim(tgfm_sig_df)[1] - aa
	n_tgfm_sig = aa +bb

	cc = sum((tgfm_null_df$sldsc_p < 0.05) & (tgfm_null_df$sldsc_z > 0.0))
	dd = dim(tgfm_null_df)[1] - cc
	n_tgfm_null = cc + dd

	fdr_vec <- c(fdr_vec, c(paste0("FDR <= ", fdr), paste0("FDR > ",0.05)))
	replication_frac_vec = c(replication_frac_vec, c(aa/n_tgfm_sig, cc/n_tgfm_null))
	totals_vec <- c(totals_vec, c(paste0(n_tgfm_sig), paste0(n_tgfm_null)))
	pval_vec <- c(pval_vec, c(".05", ".05"))


	# Now split up data frame into trait-tissue pairs identified as significant by TGFM
	tgfm_sig_df <- df[as.character(df$tgfm_sig) == "True",]
	tgfm_null_df <- df[as.character(df$tgfm_sig) == "False",]

	aa = sum((tgfm_sig_df$sldsc_p < 0.001) & (tgfm_sig_df$sldsc_z > 0.0))
	bb = dim(tgfm_sig_df)[1] - aa
	n_tgfm_sig = aa +bb

	cc = sum((tgfm_null_df$sldsc_p < 0.001) & (tgfm_null_df$sldsc_z > 0.0))
	dd = dim(tgfm_null_df)[1] - cc
	n_tgfm_null = cc + dd

	fdr_vec <- c(fdr_vec, c(paste0("FDR <= ", fdr), paste0("FDR > ",fdr)))
	replication_frac_vec = c(replication_frac_vec, c(aa/n_tgfm_sig, cc/n_tgfm_null))
	totals_vec <- c(totals_vec, c(paste0(n_tgfm_sig), paste0(n_tgfm_null)))
	pval_vec <- c(pval_vec, c(".001", ".001"))

	plot_df <- data.frame(tgfm_fdr=fdr_vec, replication_rate=replication_frac_vec, total=totals_vec, pval=factor(pval_vec))

	red_color =brewer.pal(n = 9, name = "Reds")[6]
	red_color2 =brewer.pal(n = 9, name = "Reds")[4]

	print(plot_df)


	pp <- ggplot(data=plot_df, aes(x=pval, y=replication_rate, fill=tgfm_fdr)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_text(aes(label=total, y=c(0.033,0.033, .033, .033)),vjust=0.0, position = position_dodge(0.9), color="white") +
  		labs(x="SLDSC p-value threshold", y="Replication\n rate in SLDSC") +
  		figure_theme() +
  		scale_fill_manual(values=c(red_color,red_color2)) + 
  		theme(legend.position="bottom") +
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) 
  	return(pp)
}


make_ldl_silver_standard_stratefied_gene_pip_ecdf_plot <- function(ldl_enrichment_file) {
	df <- read.table(ldl_enrichment_file, header=TRUE, sep="\t")

	df$silver_standard = factor(df$silver_standard, levels=c("silver standard", "background"))


	green_color=brewer.pal(n = 9, name = "Greens")[6]


	pp <- ggplot(df, aes(x=gene_pip, colour = silver_standard,fill=silver_standard)) +
  		stat_ecdf() +
  		figure_theme() + 
  		scale_fill_manual(values=c(green_color, "grey32")) +
  		scale_colour_manual(values=c(green_color, "grey32")) +
  		theme(legend.position="top") +
  		labs(x="TGFM (Gene) PIP", colour="", fill="") 
  	return(pp)

}

make_ldl_silver_standard_stratefied_gene_pip_density_plot <- function(ldl_enrichment_file) {
	df <- read.table(ldl_enrichment_file, header=TRUE, sep="\t")

	df$silver_standard = factor(df$silver_standard, levels=c("silver standard", "background"))


	green_color=brewer.pal(n = 9, name = "Greens")[6]


	pp <- ggplot(df, aes(x=gene_pip, colour = silver_standard,fill=silver_standard)) +
  		geom_density(alpha = 0.2) +
  		figure_theme() + 
  		scale_fill_manual(values=c(green_color, "grey32")) +
  		scale_colour_manual(values=c(green_color, "grey32")) +
  		theme(legend.position="top") +
  		labs(x="TGFM (Gene) PIP", colour="", fill="") 
  	return(pp)

}

make_ldl_silver_standard_stratefied_gene_pip_histogram <- function(ldl_enrichment_file) {
	df <- read.table(ldl_enrichment_file, header=TRUE, sep="\t")

	df$silver_standard = factor(df$silver_standard, levels=c("silver standard", "background"))


	green_color=brewer.pal(n = 9, name = "Greens")[6]


	pp <- ggplot(df, aes(x=gene_pip, colour = silver_standard,fill=silver_standard)) +
  		geom_histogram(aes(y=..density..),alpha = 0.2, position='identity') +
  		figure_theme() + 
  		scale_fill_manual(values=c(green_color, "grey32")) +
  		scale_colour_manual(values=c(green_color, "grey32")) +
  		theme(legend.position="top") +
  		labs(x="TGFM (Gene) PIP", colour="", fill="") 
  	return(pp)

}

statistical_significance_of_power_difference_between_tgfm_pip_and_ctwas_pip_when_matched_for_fdr <- function(ldl_tgfm_enrichment_file, ldl_ctwas_enrichment_file, pip_threshold, pip_lb=0.01) {
	n_total_pos = 69
	pip_vec <- c()
	efdr_vec <- c()
	efdr_lb_vec <- c()
	efdr_ub_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()
	method_vec <- c()

	# TGFM
	method_name <- "TGFM (Gene)"
	df <- read.table(ldl_tgfm_enrichment_file, header=TRUE, sep="\t")
	unique_gene_pips <- sort(unique(df$gene_pip))
	unique_gene_pips = rev(unique_gene_pips[unique_gene_pips > pip_lb])

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))


		power = n_silver/n_total_pos
		power_variance = (power*(1.0-power))/(n_total_pos)
		power_se = sqrt(power_variance)
		power_vec <- c(power_vec, power)
		power_lb_vec <- c(power_lb_vec, power - (1.96*power_se))
		power_ub_vec <- c(power_ub_vec, power + (1.96*power_se))
		method_vec <- c(method_vec, method_name)
	}


	# cTWAS
	method_name <- "cTWAS"
	df <- read.table(ldl_ctwas_enrichment_file, header=TRUE, sep="\t")
	unique_gene_pips <- sort(unique(df$gene_pip))
	unique_gene_pips = rev(unique_gene_pips[unique_gene_pips > pip_lb])

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))


		power = n_silver/n_total_pos
		power_variance = (power*(1.0-power))/(n_total_pos)
		power_se = sqrt(power_variance)
		power_vec <- c(power_vec, power)
		power_lb_vec <- c(power_lb_vec, power - (1.96*power_se))
		power_ub_vec <- c(power_ub_vec, power + (1.96*power_se))
		method_vec <- c(method_vec, method_name)
	}


	df2 <- data.frame(eFDR=efdr_vec, power=power_vec, method=factor(method_vec, levels=c("TGFM (Gene)", "cTWAS")), PIP=pip_vec, eFDR_lb=efdr_lb_vec, eFDR_ub=efdr_ub_vec, power_lb=power_lb_vec, power_ub=power_ub_vec)


	print(head(df2))

}

make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_emperical_power_curve <- function(ldl_tgfm_enrichment_file, ldl_ctwas_enrichment_file, pip_lb=0.01) {
	n_total_pos = 69
	pip_vec <- c()
	efdr_vec <- c()
	efdr_lb_vec <- c()
	efdr_ub_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()
	method_vec <- c()

	# TGFM
	method_name <- "TGFM (Gene)"
	df <- read.table(ldl_tgfm_enrichment_file, header=TRUE, sep="\t")
	unique_gene_pips <- sort(unique(df$gene_pip))
	unique_gene_pips = rev(unique_gene_pips[unique_gene_pips > pip_lb])

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))


		power = n_silver/n_total_pos
		power_variance = (power*(1.0-power))/(n_total_pos)
		power_se = sqrt(power_variance)
		power_vec <- c(power_vec, power)
		power_lb_vec <- c(power_lb_vec, power - (1.96*power_se))
		power_ub_vec <- c(power_ub_vec, power + (1.96*power_se))
		method_vec <- c(method_vec, method_name)
	}


	# cTWAS
	method_name <- "cTWAS"
	df <- read.table(ldl_ctwas_enrichment_file, header=TRUE, sep="\t")
	unique_gene_pips <- sort(unique(df$gene_pip))
	unique_gene_pips = rev(unique_gene_pips[unique_gene_pips > pip_lb])

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))


		power = n_silver/n_total_pos
		power_variance = (power*(1.0-power))/(n_total_pos)
		power_se = sqrt(power_variance)
		power_vec <- c(power_vec, power)
		power_lb_vec <- c(power_lb_vec, power - (1.96*power_se))
		power_ub_vec <- c(power_ub_vec, power + (1.96*power_se))
		method_vec <- c(method_vec, method_name)
	}


	df2 <- data.frame(eFDR=efdr_vec, power=power_vec, method=factor(method_vec, levels=c("TGFM (Gene)", "cTWAS")), PIP=pip_vec, eFDR_lb=efdr_lb_vec, eFDR_ub=efdr_ub_vec, power_lb=power_lb_vec, power_ub=power_ub_vec)

	green_color = brewer.pal(n = 9, name = "Greens")[6]
	orange_color = brewer.pal(n = 9, name = "Oranges")[6]

	pp <- ggplot(df2, aes(x = power, y = eFDR, color=method)) +
	 geom_line(size=1.0) +
	 figure_theme() + 
	 scale_colour_manual(values = c(green_color, orange_color)) +
	 labs(x="Power", y="FDR",color="") +
	 theme(legend.position="bottom")

	return(pp)

}

make_ldl_silver_standard_stratefied_gene_pip_one_minus_emperical_fdr_plot <- function(ldl_enrichment_file, pip_lb=0.01, color='green', method_name="TGFM") {
	df <- read.table(ldl_enrichment_file, header=TRUE, sep="\t")

	unique_gene_pips <- sort(unique(df$gene_pip))
	unique_gene_pips = unique_gene_pips[unique_gene_pips > pip_lb]

	pip_vec <- c()
	efdr_vec <- c()
	efdr_lb_vec <- c()
	efdr_ub_vec <- c()

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))
	}

	efdr_lb_vec[efdr_lb_vec < 0.0] = 0.0
	efdr_ub_vec[efdr_lb_vec > 1.0] = 0.0

	if (color == "green") {
		color_choice=brewer.pal(n = 9, name = "Greens")[6]
	}
	if (color == "orange") {
		color_choice=brewer.pal(n = 9, name = "Oranges")[6]
	}
	df2 <- data.frame(PIP=pip_vec, proportion=100*(1-efdr_vec), proportion_ub = 100*(1.0-efdr_lb_vec), propotion_lb=100*(1.0-efdr_ub_vec))
	pp <- ggplot(df2, aes(x = PIP, y = proportion, ymin=propotion_lb, ymax=proportion_ub)) +
	 geom_line(color=color_choice, size=2.0) +
	 geom_ribbon(alpha=0.5, fill=color_choice) +
	 figure_theme() + 
	 geom_abline(slope=100,size=.3) +
	 labs(x=paste0(method_name," (Gene) PIP"), y="% in LDL gene set  ")

	 return(pp) 


}

make_ldl_silver_standard_stratefied_gene_tissue_pip_emperical_fdr_plot <- function(ldl_enrichment_file, pip_lb=0.01, color='red', method_name="TGFM",expected_fdr=FALSE) {
	df <- read.table(ldl_enrichment_file, header=TRUE, sep="\t")

	unique_gene_pips <- sort(unique(df$gene_tissue_pip))
	unique_gene_pips = unique_gene_pips[unique_gene_pips > pip_lb]

	pip_vec <- c()
	efdr_vec <- c()
	efdr_lb_vec <- c()
	efdr_ub_vec <- c()
	expected_fdr_vec <- c()

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_tissue_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))

		expected_fdr_vec <- c(expected_fdr_vec, 1.0 - mean(tmp_df$gene_tissue_pip))
	}

	efdr_lb_vec[efdr_lb_vec < 0.0] = 0.0
	efdr_ub_vec[efdr_lb_vec > 1.0] = 0.0

	if (color == "red") {
		color_choice=brewer.pal(n = 9, name = "Reds")[6]
	}
	if (color == "green") {
		color_choice=brewer.pal(n = 9, name = "Greens")[6]
	}
	if (color == "orange") {
		color_choice=brewer.pal(n = 9, name = "Oranges")[6]
	}
	df2 <- data.frame(PIP=pip_vec, efdr=efdr_vec, efdr_lb = efdr_lb_vec, efdr_ub=efdr_ub_vec, expected_fdr=expected_fdr_vec)
	pp <- ggplot(df2, aes(x = PIP, y = efdr, ymin=efdr_lb, ymax=efdr_ub)) +
	 geom_line(color=color_choice, size=2.0) +
	 geom_ribbon(alpha=0.2, fill=color_choice) +
	 figure_theme() + 
	 geom_abline(slope=-1, intercept=1, size=.3) +
	 labs(x=paste0(method_name," (Gene-Tissue) PIP"), y="FDR (Silver standard gene-tissue set)")

	 if (expected_fdr==TRUE) {
	 	pp <- pp + geom_line(data=df2, aes(x=PIP, y=expected_fdr), linetype="dashed", colour="grey30",size=1.3)

	 }


	 return(pp) 


}

get_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_data <- function(ldl_enrichment_file, pip_lb=0.01, color='green', method_name="TGFM") {
	df <- read.table(ldl_enrichment_file, header=TRUE, sep="\t")

	unique_gene_pips <- sort(unique(df$gene_pip))
	unique_gene_pips = unique_gene_pips[unique_gene_pips > pip_lb]

	pip_vec <- c()
	efdr_vec <- c()
	efdr_lb_vec <- c()
	efdr_ub_vec <- c()
	expected_fdr_vec <- c()

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))

		expected_fdr_vec <- c(expected_fdr_vec, 1.0 - mean(tmp_df$gene_pip))
	}

	efdr_lb_vec[efdr_lb_vec < 0.0] = 0.0
	efdr_ub_vec[efdr_lb_vec > 1.0] = 0.0

	df2 <- data.frame(PIP=pip_vec, efdr=efdr_vec, efdr_lb = efdr_lb_vec, efdr_ub=efdr_ub_vec,expected_fdr=expected_fdr_vec)

	return(df2)


}

statistical_significance_of_fdr_difference_between_tgfm_pip_and_ctwas_pip_calibration <- function(tgfm_ldl_enrichment_file, ctwas_ldl_enrichment_file, pip_threshold) {
	df_tgfm <- read.table(tgfm_ldl_enrichment_file, header=TRUE, sep="\t")
	df_ctwas <- read.table(ctwas_ldl_enrichment_file, header=TRUE, sep="\t")


	# TGFM calc
	tgfm_tmp_df <- df_tgfm[df_tgfm$gene_pip >= pip_threshold,]
	tgfm_n_silver = sum(as.character(tgfm_tmp_df$silver_standard) == "silver standard")
	tgfm_n_background = sum(as.character(tgfm_tmp_df$silver_standard) == "background")

	# cTWAS calc
	ctwas_tmp_df <- df_ctwas[df_ctwas$gene_pip >= pip_threshold,]
	ctwas_n_silver = sum(as.character(ctwas_tmp_df$silver_standard) == "silver standard")
	ctwas_n_background = sum(as.character(ctwas_tmp_df$silver_standard) == "background")

	print(tgfm_n_silver)
	print(tgfm_n_background)

	print(ctwas_n_silver)
	print(ctwas_n_background)


	data <- matrix(c(tgfm_n_silver, tgfm_n_background, ctwas_n_silver, ctwas_n_background), nrow = 2)

	# Perform Fisher's exact test
	result <- fisher.test(data)
	print(result)

}


make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_plot <- function(ldl_enrichment_file, pip_lb=0.01, color='green', method_name="TGFM", expected_fdr=FALSE) {
	df <- read.table(ldl_enrichment_file, header=TRUE, sep="\t")

	unique_gene_pips <- sort(unique(df$gene_pip))
	unique_gene_pips = unique_gene_pips[unique_gene_pips > pip_lb]

	pip_vec <- c()
	efdr_vec <- c()
	efdr_lb_vec <- c()
	efdr_ub_vec <- c()
	expected_fdr_vec <- c()

	for (pip_iter in 1:length(unique_gene_pips)) {
		pip_val = unique_gene_pips[pip_iter]
		tmp_df <- df[df$gene_pip >= pip_val,]

		n_silver = sum(as.character(tmp_df$silver_standard) == "silver standard")
		n_background = sum(as.character(tmp_df$silver_standard) == "background")

		efdr = n_background/(n_silver + n_background)

		efdr_variance = (efdr*(1.0-efdr))/(n_silver + n_background)
		efdr_se = sqrt(efdr_variance)

		pip_vec <- c(pip_vec, pip_val)
		efdr_vec <- c(efdr_vec, efdr)

		efdr_lb_vec <- c(efdr_lb_vec, efdr - (1.96*efdr_se))
		efdr_ub_vec <- c(efdr_ub_vec, efdr + (1.96*efdr_se))


		expected_fdr_vec <- c(expected_fdr_vec, 1.0 - mean(tmp_df$gene_pip))
	}

	efdr_lb_vec[efdr_lb_vec < 0.0] = 0.0
	efdr_ub_vec[efdr_lb_vec > 1.0] = 0.0

	if (color == "green") {
		color_choice=brewer.pal(n = 9, name = "Greens")[6]
	}
	if (color == "orange") {
		color_choice=brewer.pal(n = 9, name = "Oranges")[6]
	}
	df2 <- data.frame(PIP=pip_vec, efdr=efdr_vec, efdr_lb = efdr_lb_vec, efdr_ub=efdr_ub_vec, expected_fdr=expected_fdr_vec)
	pp <- ggplot(df2, aes(x = PIP, y = efdr, ymin=efdr_lb, ymax=efdr_ub)) +
	 geom_line(color=color_choice, size=2.0) +
	 geom_ribbon(alpha=0.2, fill=color_choice) +
	 figure_theme() + 
	 geom_abline(slope=-1, intercept=1, size=.3) +
	 labs(x=paste0(method_name," (Gene) PIP"), y="FDR")

	 if (expected_fdr==TRUE) {
	 	pp <- pp + geom_line(data=df2, aes(x=PIP, y=expected_fdr), linetype="dashed", colour="grey30",size=1.3) 

	 }

	 return(pp) 


}



##################################
# Extract command line arguments
##################################

independent_trait_names_file <- args[1]
tgfm_results_dir <- args[2]
tgfm_organized_results_dir <- args[3]
gtex_tissue_colors_file <- args[4]
iterative_tgfm_prior_dir <- args[5]
pops_enrichment_dir <- args[6]
non_disease_specific_gene_set_enrichment_dir <- args[7]
visualize_tgfm_dir <- args[8]
tissue_tissue_correlation_file <- args[9]
gtex_pseudotissue_file <- args[10]
chromatin_cell_type_group_ldsc_dir <- args[11]
ldl_silver_standard_gene_set_enrichment_dir <- args[12]

##########################################################
# Load in data
##########################################################
# Independent traits used in https://www.nature.com/articles/s41588-020-00735-5
independent_traits <- c("body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled")

# Load in pseudotissue names and sample sizes
pseudotissue_df <- read.table(gtex_pseudotissue_file, header=TRUE, sep="\t")
pseudotissue_names <- pseudotissue_df$pseudotissue_name
pseudotissue_ss <- pseudotissue_df$sample_size

# Load in gtex tissue colors 
gtex_colors_df <- read.table(gtex_tissue_colors_file, header=TRUE, sep="\t")
gtex_colors_df$tissue_site_detail_id = as.character(gtex_colors_df$tissue_site_detail_id)
gtex_colors_df$tissue_site_detail_id[23] = "Cells_Cultured_fibroblasts"

# Extract trait names
trait_df <- read.table(independent_trait_names_file, header=TRUE, sep="\t")
trait_names <- as.character(trait_df$study_name)
trait_names_readable <- as.character(trait_df$study_name_readable)



gene_type="component_gene"
# Extract tissue names
iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_names[1], "_", gene_type,"_", "susie_pmces_uniform", "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
aa = read.table(iterative_prior_file, header=TRUE,sep="\t")
tissue_names = as.character(aa$element_name[2:length(aa$element_name)])



##################################################
# Make tissue-tissue correlation heatmap
##################################################
if (FALSE) {
output_file <- paste0(visualize_tgfm_dir, "tissue_tissue_correlation_heatmap.pdf")
tt_corr_heatmap = make_tissue_tissue_correlation_heatmap(tissue_tissue_correlation_file)
ggsave(tt_corr_heatmap, file=output_file, width=7.2, height=7.5, units="in")


##################################################
# Make scatterplot showing relationship between sample size and average correlation
##################################################
output_file <- paste0(visualize_tgfm_dir, "tissue_tissue_correlation_sample_size_scatter.pdf")
tt_corr_ss_scatter = make_tissue_tissue_correlation_sample_size_scatterplot(tissue_tissue_correlation_file, pseudotissue_names, pseudotissue_ss)
ggsave(tt_corr_ss_scatter, file=output_file, width=7.2, height=5.5, units="in")
}

##################################################
# Simply count up number of gene-tissue pairs, genes, and variants above some PIP threshold
##################################################
if (FALSE) {
pip_threshold = "0.5"
count_up_number_of_genetic_elements(trait_names, tgfm_organized_results_dir,pip_threshold)
pip_threshold = "0.9"
count_up_number_of_genetic_elements(trait_names, tgfm_organized_results_dir,pip_threshold)
}

##################################################
# standard error barplot showing non-disease-specific gene set enrichments
##################################################
if (FALSE) {
enrichment_barplot <- make_non_disease_specific_gene_set_enrichment_barplot_at_single_pip_thresh(non_disease_specific_gene_set_enrichment_dir)
output_file <- paste0(visualize_tgfm_dir, "non_disease_specific_gene_set_enrichments_standard_error_barplot.pdf")
ggsave(enrichment_barplot, file=output_file, width=7.2, height=4.6, units="in")

# X threshold
x_thresh_enrichment_barplot <- make_non_disease_specific_gene_set_enrichment_barplot_cross_pip_thresholds(non_disease_specific_gene_set_enrichment_dir)
output_file <- paste0(visualize_tgfm_dir, "non_disease_specific_gene_set_enrichments_cross_pip_thresholds_standard_error_barplot.pdf")
ggsave(x_thresh_enrichment_barplot, file=output_file, width=7.2, height=4.6, units="in")


# Joint plot
joint_enrichment_barplot <- plot_grid(enrichment_barplot, x_thresh_enrichment_barplot+theme(legend.position="top"), ncol=2, labels=c("a","b"))
output_file <- paste0(visualize_tgfm_dir, "non_disease_specific_gene_set_enrichments_joint_standard_error_barplot.pdf")
ggsave(joint_enrichment_barplot, file=output_file, width=7.2, height=4.0, units="in")


##################################################
# Supp data table showing standard error barplot showing non-disease-specific gene set enrichments
##################################################
supp_table_df <- make_non_disease_specific_gene_set_enrichment_suppData_cross_pip_thresh(non_disease_specific_gene_set_enrichment_dir)
supp_table_file = paste0(visualize_tgfm_dir, "suppTable_non_disease_gene_set_enrichment.txt")
write.table(supp_table_df, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
print(supp_table_file)
}




##################################################
# Violin plot showing bootstrapped prior distributions 
##################################################
if (FALSE) {
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	# Sampler approach
	gene_type="component_gene"
	method_version="susie_pmces_uniform"
	iterative_prior_file <- paste0(iterative_tgfm_prior_dir, "tgfm_results_", trait_name, "_", gene_type, "_", method_version, "_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt")
	iterative_sampler_violin_plot <- make_violin_plot_showing_distribution_of_iterative_prior_probability_of_each_tissue(iterative_prior_file, method_version, trait_name_readable)

	# Iterative version
	output_file <- paste0(visualize_tgfm_dir, "tissue_violinplot_of_distribution_prior_probabilities_", gene_type, "_", trait_name_readable,".pdf")
	ggsave(iterative_sampler_violin_plot, file=output_file, width=10.2, height=4.6, units="in")
}
}

##################################################
# Make file summarizing prior significance of each trait-tissue pair
##################################################
method_version="susie_pmces_uniform"
gene_type <- "component_gene"
trait_tissue_prior_significance_file <- paste0(visualize_tgfm_dir, gene_type, "_trait_tissue_prior_bonferronni_corrected_significance.txt")
if (FALSE) {
generate_file_containing_bonf_significance_of_each_trait_tissue_pair_based_on_iterative_prior(trait_names, iterative_tgfm_prior_dir, method_version,trait_tissue_prior_significance_file, gene_type)
}




##########################################################
# Make TGFM PIP density across genetic element classes
##########################################################
if (FALSE) {
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_tgfm_pip_density_plot_stratified_by_element_class.pdf")
density_hist <- make_pip_density_histogram_across_element_classes(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
ggsave(density_hist, file=output_file_pdf, width=7.2, height=3.7, units="in", dpi=400)


##########################################################
# Swarm-plot showing gene-tissue PIPs colored by tissue group for each trait
##########################################################
gene_type <- "component_gene"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
beeswarm_plot <- make_swarm_plot_showing_gene_tissue_pips_colored_by_tissue_group_for_each_trait(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, gene_type)
output_file <- paste0(visualize_tgfm_dir, "beeswarm_gene_tissue_pip_colored_by_tissue_group_", method_version, "_", gene_type, ".pdf")
ggsave(beeswarm_plot, file=output_file, width=15.2, height=4.5, units="in")

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
med_prob_se_barplot2 <- make_expected_fraction_of_genetic_elements_from_gene_expression_se_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
ggsave(med_prob_se_barplot2, file=output_file, width=7.2, height=4.2, units="in")

##########################################################
# Joint fraction mediated plot
##########################################################
joint_fraction_med_plot <- plot_grid(plot_grid(NULL,med_prob_se_barplot, ncol=2,rel_widths=c(.02,1)),plot_grid(NULL,med_prob_se_barplot2,ncol=2,rel_widths=c(.02,1.0)), ncol=1, labels=c("a","b"), rel_heights=c(1.3,1.0))
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_joint_fraction_mediated_plot.pdf")
ggsave(joint_fraction_med_plot, file=output_file, width=7.2, height=6.2, units="in")
}


##########################################################
# Bar plot showing expected number of causal gene-tissue pairs in each tissue
##########################################################
if (FALSE) {
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	trait_name_readable <- trait_names_readable[trait_iter]

	gene_type <- "component_gene"
	# Sampler approach
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	tissue_bar_plot <- make_bar_plot_showing_expected_number_of_causal_gene_tissue_pairs_for_single_trait(trait_name, trait_name_readable, method_version, tgfm_organized_results_dir, tissue_names, gene_type)
	output_file <- paste0(visualize_tgfm_dir, "tissue_barplot_of_expected_number_of_gene_tissue_pairs_", trait_name_readable, "_", gene_type, "_", method_version,".pdf")
	ggsave(tissue_bar_plot, file=output_file, width=7.2, height=3.0, units="in")
}
}



if (FALSE) {
##########################################################
# Barplot showing number high pip gene-tissue pairs between methods
##########################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
pip_threshold <- "0.5"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
n_genes_barplot_5 <- make_number_of_high_pip_gene_tissue_pairs_method_comparison_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits)
pip_threshold <- "0.9"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
n_genes_barplot_9 <- make_number_of_high_pip_gene_tissue_pairs_method_comparison_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits)

n_genes_joint_plot <- plot_grid(n_genes_barplot_5, n_genes_barplot_9,ncol=1)
output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_method_comparison_barplot.pdf")
ggsave(n_genes_joint_plot, file=output_file, width=7.2, height=5.7, units="in")


##########################################################
# Barplot showing number high pip genes vs gene-tissue pairs
##########################################################
pip_threshold <- "0.5"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
n_genes_barplot_5 <- make_number_of_high_pip_genes_vs_gene_tissues_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits)
pip_threshold <- "0.9"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
n_genes_barplot_9 <- make_number_of_high_pip_genes_vs_gene_tissues_barplot(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, pip_threshold, independent_traits)
# Stack together
n_genes_joint_plot <- plot_grid(n_genes_barplot_5, n_genes_barplot_9,ncol=1)

output_file <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_genes_vs_gene_tissues_barplot.pdf")
ggsave(n_genes_joint_plot, file=output_file, width=7.2, height=5.7, units="in")
}


##########################################################
# Heatmap showing expected number of causal genes in each tissue-trait pair
##########################################################
if (FALSE) {
pip_threshs <- c(.01, .1, .3, .5)
gene_type="component_gene"
for (pip_iter in 1:length(pip_threshs)) {
	pip_thresh <- pip_threshs[pip_iter]
	print(pip_thresh)
	method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
	heatmap <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_thresh, trait_tissue_prior_significance_file, gene_type)
	output_file <- paste0(visualize_tgfm_dir, "expected_num_causal_genes_", pip_thresh, "_", method_version, "_", gene_type,"_heatmap.pdf")
	ggsave(heatmap, file=output_file, width=7.2, height=8.0, units="in")
}
}


##########################################################
# Supplementary data for heatmap of expected number of causal genes in each tissue-trait pir
##########################################################
if (FALSE) {
pip_thresh=.5
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
supp_table_df <- get_heatmap_data_showing_expected_number_of_causal_gene_tissue_pairs_cross_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_thresh, trait_tissue_prior_significance_file)
supp_table_file = paste0(visualize_tgfm_dir, "suppTable_figure4a_numerical.txt")
write.table(supp_table_df, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
}

if (FALSE) {
##########################################################
# Make Figure 3
##########################################################
# Get ordered traits according to number of hits identified by gene-tissue pairs
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
ordered_traits <- make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_trait_order_according_to_gene_tissue_pairs(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)

# Make heatmap-barplot showing expected number of causal variants
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
#output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_variants_heatmap_barplot_v2.pdf")
variant_heatmap_barplot <- make_number_of_high_pip_variants_heatmap_barplot_v2(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
#ggsave(variant_heatmap_barplot, file=output_file_pdf, width=7.2, height=3.7, units="in", dpi=400)
variant_heatmap_barplot_ord <- make_number_of_high_pip_variants_heatmap_barplot_v2(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, preordered=TRUE, ordered_traits=ordered_traits)


# Make heatmap-barplot showing expected number of causal genes
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_genes_heatmap_barplot_v2.pdf")
gene_heatmap_barplot <- make_number_of_high_pip_genes_heatmap_barplot_v2(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
#ggsave(gene_heatmap_barplot, file=output_file_pdf, width=7.2, height=3.7, units="in", dpi=400)
gene_heatmap_barplot_ord <- make_number_of_high_pip_genes_heatmap_barplot_v2(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, preordered=TRUE, ordered_traits=ordered_traits)


# Make heatmap-barplot showing expected number of causal gene-tissue pairs
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
#output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2.pdf")
gene_tissue_heatmap_barplot <- make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits)
#ggsave(gene_tissue_heatmap_barplot, file=output_file_pdf, width=7.2, height=3.7, units="in", dpi=400)
gene_tissue_heatmap_barplot_ord <- make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, independent_traits, preordered=TRUE, ordered_traits=ordered_traits)

# Make joint heatmap barplot plot for different genetic element categories
red_color=brewer.pal(n = 9, name = "Blues")[8]
variant_heatmap_barplot2 <- plot_grid(NULL, variant_heatmap_barplot_ord + annotate("text",  x=Inf, y = Inf, label = "Variant", vjust=1, hjust=1, color=red_color), rel_widths=c(.02,1),ncol=2)
green_color=brewer.pal(n = 9, name = "Greens")[8]
gene_heatmap_barplot2 <- plot_grid(NULL, gene_heatmap_barplot_ord + annotate("text",  x=Inf, y = Inf, label = "Gene", vjust=1, hjust=1,color=green_color)+ theme(axis.text.x = element_blank()), rel_widths=c(.02,1),ncol=2)
blue_color=brewer.pal(n = 9, name = "Reds")[8]
gene_tissue_heatmap_barplot2 <- plot_grid(NULL, gene_tissue_heatmap_barplot_ord+ annotate("text",  x=Inf, y = Inf, label = "Gene-Tissue", vjust=1, hjust=1,color=blue_color)+ theme(axis.text.x = element_blank()), rel_widths=c(.02,1),ncol=2)
joint_heatmap_barplot <- plot_grid(NULL,gene_tissue_heatmap_barplot2, gene_heatmap_barplot2, variant_heatmap_barplot2, ncol=1, labels=c("","a", "b", "c"), rel_heights=c(.03, 1, 1,1.65))
output_file_pdf <- paste0(visualize_tgfm_dir, "figure3.pdf")
ggsave(joint_heatmap_barplot, file=output_file_pdf, width=7.2, height=6.4, units="in", dpi=400)

output_file_pdf <- paste0(visualize_tgfm_dir, "figure3_for_poster.pdf")
gene_tissue_heatmap_barplot2 <- plot_grid(NULL, gene_tissue_heatmap_barplot_ord+ annotate("text",  x=Inf, y = Inf, label = "Gene-Tissue", vjust=1, hjust=1,color=blue_color), rel_widths=c(.02,1),ncol=2)
gene_heatmap_barplot2 <- plot_grid(NULL, gene_heatmap_barplot_ord + annotate("text",  x=Inf, y = Inf, label = "Gene", vjust=1, hjust=1,color=green_color), rel_widths=c(.02,1),ncol=2)
variant_heatmap_barplot2 <- plot_grid(NULL, variant_heatmap_barplot_ord + annotate("text",  x=Inf, y = Inf, label = "Variant", vjust=1, hjust=1, color=red_color), rel_widths=c(.02,1),ncol=2)
joint_heatmap_barplot_poster <- plot_grid(gene_tissue_heatmap_barplot2, gene_heatmap_barplot2, variant_heatmap_barplot2, ncol=3)

ggsave(joint_heatmap_barplot_poster, file=output_file_pdf, width=21.2, height=2.7, units="in", dpi=400)
}

if (FALSE) {
#######################################################
# Supp figure for figure 3 across all traits
#######################################################
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
ordered_traits <- rev(make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_trait_order_according_to_gene_tissue_pairs(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, trait_names, gene_type=gene_type))
# Make heatmap-barplot showing expected number of causal variants

# Variants
variant_heatmap_t_barplot <- make_number_of_high_pip_variants_heatmap_barplot_v2_transpose(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, trait_names,gene_type=gene_type, preordered=TRUE,ordered_traits=ordered_traits) +theme(axis.title.y=element_blank(), axis.text.y=element_blank())+labs(x="No. fine-mapped\nvariants")
# Genes
gene_heatmap_t_barplot <- make_number_of_high_pip_genes_heatmap_barplot_v2_transpose(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, trait_names, gene_type=gene_type,preordered=TRUE,ordered_traits=ordered_traits)+theme(axis.title.y=element_blank(), axis.text.y=element_blank())+labs(x="No. fine-mapped\ngenes")
# Gene-tissue pairs
gene_tissue_heatmap_t_barplot <- make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2_transpose(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, trait_names, gene_type=gene_type, preordered=TRUE,ordered_traits=ordered_traits) + labs(x="No. fine-mapped\ngene-tissue pairs")


output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_", gene_type, "_number_of_high_pip_genetic_elements_heatmap_barplot_all_traits.pdf")
joint_plot <- plot_grid(gene_tissue_heatmap_t_barplot + theme(legend.position="bottom"), gene_heatmap_t_barplot + theme(legend.position="bottom"),variant_heatmap_t_barplot + theme(legend.position="bottom"), ncol=3,rel_widths=c(1.11,.6,.6))

ggsave(joint_plot, file=output_file_pdf, width=7.2, height=7.2, units="in", dpi=400)


#######################################################
# Supp figure for figure 3 (one for variant, gene, and gene-tissue)
#######################################################
# Make heatmap-barplot showing expected number of causal variants
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"

# Variants
variant_heatmap_t_barplot <- make_number_of_high_pip_variants_heatmap_barplot_v2_transpose(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, trait_names)
output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_variants_heatmap_barplot_all_traits.pdf")
ggsave(variant_heatmap_t_barplot, file=output_file_pdf, width=7.2, height=6.0, units="in", dpi=400)

# Genes
gene_heatmap_t_barplot <- make_number_of_high_pip_genes_heatmap_barplot_v2_transpose(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, trait_names)
output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_genes_heatmap_barplot_all_traits.pdf")
ggsave(gene_heatmap_t_barplot, file=output_file_pdf, width=7.2, height=6.0, units="in", dpi=400)

# Gene-tissue pairs
gene_tissue_heatmap_t_barplot <- make_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_v2_transpose(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, trait_names)
output_file_pdf <- paste0(visualize_tgfm_dir, "tgfm_", method_version, "_number_of_high_pip_gene_tissue_pairs_heatmap_barplot_all_traits.pdf")
ggsave(gene_tissue_heatmap_t_barplot, file=output_file_pdf, width=7.2, height=6.0, units="in", dpi=400)
}



##########################################################
# Make heatmap (with significance stars) showing expected number of causal genes in each tissue-trait pair for selected traits
##########################################################
if (FALSE) {
pip_thresh <- "0.5"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
selected_traits <- c("disease_AID_ALL", "biochemistry_VitaminD", "body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled")
heatmap <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_thresh, selected_traits,trait_tissue_prior_significance_file, significance_bool=TRUE, only_fdr_05=TRUE)
output_file <- paste0(visualize_tgfm_dir, "expected_num_causal_genes_", pip_thresh, "_", method_version,"_selected_traits_heatmap_with_significance.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")
}

##########################################################
# Make heatmap showing expected number of causal genes in each tissue-trait pair for selected traits
##########################################################
if (FALSE) {
pip_thresh <- "0.5"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
selected_traits <- c("disease_AID_ALL", "biochemistry_VitaminD", "body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled")
heatmap <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_thresh, selected_traits,trait_tissue_prior_significance_file)
output_file <- paste0(visualize_tgfm_dir, "expected_num_causal_genes_", pip_thresh, "_", method_version,"_selected_traits_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")
}

##########################################################
# Make POPS enrichment standard error barplot
##########################################################
if (FALSE) {
# Load in summary data
pops_summary_df <- read.table(paste0(pops_enrichment_dir, "cross_traits_pops_tgfm_enrichment_summary.txt"), header=TRUE)
# average and standard error of mean of pops-score binned by TGFM PIP
output_file <- paste0(visualize_tgfm_dir, "mean_se_barplot_pops_score_binned_by_tgfm_pip_x_trait.pdf")
barplot <- mean_se_barplot_of_pops_score_binned_by_tgfm_pip(pops_summary_df, independent_traits)
ggsave(barplot, file=output_file, width=7.2, height=3.7, units="in")
# Make PoPS supp data file
pops_summary_df <- read.table(paste0(pops_enrichment_dir, "cross_traits_pops_tgfm_enrichment_summary.txt"), header=TRUE)
supp_table_df <- mean_se_barplot_of_pops_score_binned_by_tgfm_pip_supp_table(pops_summary_df, independent_traits)
supp_table_file = paste0(visualize_tgfm_dir, "suppTable_figure4c_numerical.txt")
write.table(supp_table_df, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
# Load in summary data
pops_summary_df <- read.table(paste0(pops_enrichment_dir, "cross_traits_pops_tgfm_enrichment_summary.txt"), header=TRUE)
# average and standard error of mean of pops-score binned by TGFM PIP
output_file <- paste0(visualize_tgfm_dir, "mean_se_barplot_pops_score_binned_by_tgfm_pip_x_trait.pdf")
barplot <- mean_se_barplot_of_pops_score_binned_by_tgfm_pip(pops_summary_df, independent_traits)
ggsave(barplot, file=output_file, width=7.2, height=3.7, units="in")
}

##########################################################
# Make trait-tissue pair chromatin overlap bar plot
##########################################################
if (FALSE) {
fdr <- "0.05"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot <- make_trait_tissue_chromatin_overlap_barplot(trait_tissue_chromatin_overlap_barplot, independent_traits=TRUE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_fdr_", fdr, ".pdf")
ggsave(barplot, file=output_file, width=7.2, height=4.7, units="in")
}

if (FALSE) {
##########################################################
# Make LDL silver standard stratified TGFM gene pip histogram
##########################################################
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
histo <- make_ldl_silver_standard_stratefied_gene_pip_histogram(ldl_enrichment_file)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_gene_set_tgfm_enrichment_histogram_plot.pdf")
ggsave(histo, file=output_file, width=7.2, height=4.7, units="in")
print(output_file)

##########################################################
# Make LDL silver standard stratified TGFM gene pip density
##########################################################
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
histo <- make_ldl_silver_standard_stratefied_gene_pip_density_plot(ldl_enrichment_file)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_gene_set_tgfm_enrichment_density_plot.pdf")
ggsave(histo, file=output_file, width=7.2, height=4.7, units="in")
print(output_file)

##########################################################
# Make LDL silver standard stratified TGFM gene pip eCDF
##########################################################
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
ecdf <- make_ldl_silver_standard_stratefied_gene_pip_ecdf_plot(ldl_enrichment_file)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_gene_set_tgfm_enrichment_ecdf.pdf")
ggsave(ecdf, file=output_file, width=7.2, height=4.7, units="in")
print(output_file)
}

##########################################################
# Make LDL silver standard emperical FDR plot
##########################################################
if (FALSE) {
# TGFM
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
tgfm_eFDR_plot <- make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_plot(ldl_enrichment_file)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_tgfm_gene_set_tgfm_empirical_FDR.pdf")
ggsave(tgfm_eFDR_plot, file=output_file, width=7.2, height=4.7, units="in")

# cTWAS
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_ctwas_enrichment_summary.txt")
ctwas_eFDR_plot <- make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_plot(ldl_enrichment_file,color="orange", method="cTWAS")
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_ctwas_gene_set_tgfm_empirical_FDR.pdf")
ggsave(ctwas_eFDR_plot, file=output_file, width=7.2, height=4.7, units="in")

# Joint plot 
joint_eFDR_plot <- plot_grid(tgfm_eFDR_plot, ctwas_eFDR_plot, ncol=2)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_gene_set_tgfm_empirical_FDR.pdf")
ggsave(joint_eFDR_plot, file=output_file, width=7.2, height=3.5, units="in")

# TGFM gene-tissue
ldl_gt_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_gene_tissue_enrichment_summary.txt")
tgfm_gt_eFDR_plot <- make_ldl_silver_standard_stratefied_gene_tissue_pip_emperical_fdr_plot(ldl_gt_enrichment_file)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_tgfm_gene_tissue_set_tgfm_empirical_FDR.pdf")
ggsave(tgfm_gt_eFDR_plot, file=output_file, width=7.2, height=4.7, units="in")

# Save raw data to supplementary table
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
supp_table_df <- get_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_data(ldl_enrichment_file)
supp_table_file = paste0(visualize_tgfm_dir, "suppTable_figure4d_numerical.txt")
write.table(supp_table_df, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
}

##########################################################
# Make LDL silver standard emperical FDR plot (show expected PIP)
##########################################################
# TGFM
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
tgfm_eFDR_plot <- make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_plot(ldl_enrichment_file, expected_fdr=TRUE)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_tgfm_gene_set_tgfm_empirical_FDR_w_expected_fdr.pdf")
ggsave(tgfm_eFDR_plot, file=output_file, width=7.2, height=4.7, units="in")

# cTWAS
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_ctwas_enrichment_summary.txt")
ctwas_eFDR_plot <- make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_plot(ldl_enrichment_file,color="orange", method="cTWAS", expected_fdr=TRUE)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_ctwas_gene_set_tgfm_empirical_FDR_w_expected_fdr.pdf")
ggsave(ctwas_eFDR_plot, file=output_file, width=7.2, height=4.7, units="in")

# Joint plot 
joint_eFDR_plot <- plot_grid(tgfm_eFDR_plot, ctwas_eFDR_plot, ncol=2, labels=c("a","b"))
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_gene_set_tgfm_empirical_FDR_w_expected_fdr.pdf")
ggsave(joint_eFDR_plot, file=output_file, width=7.2, height=3.5, units="in")


##########################################################
# Make LDL silver standard emperical FDR-empirical power curv
##########################################################
ldl_tgfm_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
ldl_ctwas_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_ctwas_enrichment_summary.txt")
eFDR_ePower_curve <- make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_emperical_power_curve(ldl_tgfm_enrichment_file, ldl_ctwas_enrichment_file)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_gene_set_empirical_FDR_empirical_power_curve.pdf")
ggsave(eFDR_ePower_curve, file=output_file, width=7.2, height=4.1, units="in")

# Second joint plot
joint_eFDR_plot2 <- plot_grid(joint_eFDR_plot, eFDR_ePower_curve, ncol=1, labels=c("", "c"))
output_file <- paste0(visualize_tgfm_dir, "tgfm_gene_ctwas_comparison_with_silver_standard_gene_set_supp_fig.pdf")
ggsave(joint_eFDR_plot2, file=output_file, width=7.2, height=5.8, units="in")

# TGFM gene-tissue
ldl_gt_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_gene_tissue_enrichment_summary.txt")
tgfm_gt_eFDR_plot <- make_ldl_silver_standard_stratefied_gene_tissue_pip_emperical_fdr_plot(ldl_gt_enrichment_file, expected_fdr=TRUE)
output_file <- paste0(visualize_tgfm_dir, "ldl_silver_standard_tgfm_gene_tissue_set_tgfm_empirical_FDR_w_expected_fdr.pdf")
ggsave(tgfm_gt_eFDR_plot, file=output_file, width=7.2, height=4.7, units="in")

##########################################################
# Statistical signficance of FDR difference between TGFM PIP and ctwas pip
##########################################################
tgfm_ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
ctwas_ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_ctwas_enrichment_summary.txt")
if (FALSE) {
statistical_significance_of_fdr_difference_between_tgfm_pip_and_ctwas_pip_calibration(tgfm_ldl_enrichment_file, ctwas_ldl_enrichment_file, 0.5)
}

##########################################################
# Statistical signficance of power difference between TGFM and ctwas when matched for fdr
##########################################################
tgfm_ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
ctwas_ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_ctwas_enrichment_summary.txt")
if (FALSE) {
statistical_significance_of_power_difference_between_tgfm_pip_and_ctwas_pip_when_matched_for_fdr(tgfm_ldl_enrichment_file, ctwas_ldl_enrichment_file, 0.5)
}





##########################################################
# Make trait-tissue pair chromatin overlap bar plot for a range of thresholds
##########################################################
if (FALSE) {
fdr <- "0.05"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot_05 <- make_trait_tissue_chromatin_overlap_barplot_multiple_thresh(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=TRUE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_multiple_thresh.pdf")
ggsave(barplot_05, file=output_file, width=7.2, height=4.7, units="in")
}

if (FALSE) {
##########################################################
# Make trait-tissue pair chromatin overlap histogram plot
##########################################################
fdr <- "0.05"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot_05 <- make_trait_tissue_chromatin_overlap_histogram(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_histogram_", fdr, ".pdf")
ggsave(barplot_05, file=output_file, width=7.2, height=4.7, units="in")

##########################################################
# Make trait-tissue pair chromatin overlap densityp plot
##########################################################
fdr <- "0.05"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot_05 <- make_trait_tissue_chromatin_overlap_density(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_density_", fdr, ".pdf")
ggsave(barplot_05, file=output_file, width=7.2, height=4.7, units="in")
}

##########################################################
# Make trait-tissue pair chromatin overlap cdf plot
##########################################################
fdr <- "0.05"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot_05 <- make_trait_tissue_chromatin_overlap_cdf(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_cdf_", fdr, ".pdf")
ggsave(barplot_05, file=output_file, width=7.2, height=4.7, units="in")

fdr <- "0.2"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot_2 <- make_trait_tissue_chromatin_overlap_cdf(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_cdf_", fdr, ".pdf")
ggsave(barplot_2, file=output_file, width=7.2, height=4.7, units="in")

fdr <- "0.5"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
print(trait_tissue_chromatin_overlap_barplot)
barplot_2 <- make_trait_tissue_chromatin_overlap_cdf(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_cdf_", fdr, ".pdf")
ggsave(barplot_2, file=output_file, width=7.2, height=4.7, units="in")

fdr <- "0.75"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot_2 <- make_trait_tissue_chromatin_overlap_cdf(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_cdf_", fdr, ".pdf")
ggsave(barplot_2, file=output_file, width=7.2, height=4.7, units="in")

fdr <- "0.9"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
barplot_2 <- make_trait_tissue_chromatin_overlap_cdf(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
output_file <- paste0(visualize_tgfm_dir, "trait_tissue_pair_chromatin_overlap_tgfm_cdf_", fdr, ".pdf")
ggsave(barplot_2, file=output_file, width=7.2, height=4.7, units="in")

# Get raw data for supp table
fdr <- "0.05"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
supp_table_df <- get_trait_tissue_chromatin_overlap_cdf_raw_data(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE)
supp_table_file = paste0(visualize_tgfm_dir, "suppTable_figure4b_numerical.txt")
write.table(supp_table_df, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)

##########################################################
# Make Figure 4 
##########################################################
red_color =brewer.pal(n = 9, name = "Reds")[7]
# FIG 4A
pip_thresh <- "0.5"
method_version="susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
selected_traits <- c("disease_AID_ALL", "biochemistry_VitaminD", "body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled")
fig_4a <- make_heatmap_showing_expected_number_of_causal_gene_tissue_pairs_cross_selected_traits(trait_names, trait_names_readable, method_version, tgfm_organized_results_dir, tissue_names, pip_thresh, selected_traits,trait_tissue_prior_significance_file, significance_bool=TRUE, only_fdr_05=TRUE)+ 
          theme(legend.position="top") + theme(plot.margin = margin(5.5, 5.5, 0.0, 0.0, "points")) + scale_fill_gradient(low = "white", high = red_color, guide = guide_colorbar(title.position = "top")) +
          theme(legend.margin = margin(t = 0, r = 65, b = 0, l = 0)) +
          theme(axis.text.x = element_text(angle = 45, hjust=1)) +
          labs(x = NULL, y=NULL)
# FIG 4B
fdr <- "0.05"
trait_tissue_chromatin_overlap_barplot <- paste0(chromatin_cell_type_group_ldsc_dir, "tgfm_sldsc_chromatin_overlap_summary_", fdr, ".txt")
fig_4b <- make_trait_tissue_chromatin_overlap_cdf(trait_tissue_chromatin_overlap_barplot, fdr, independent_traits=FALSE) +
          theme(legend.position="top") + guides(fill = guide_legend(title.position = "top")) + labs(fill="TGFM tissue-specific prior", colour="TGFM tissue-specific prior") +
          theme(plot.margin = margin(5.5, 5.5, 0.0, 5.5, "points"))

fig_4b <- plot_grid(fig_4b, NULL, ncol=1,rel_heights=c(1,.16))


# FIG 4c
pops_summary_df <- read.table(paste0(pops_enrichment_dir, "cross_traits_pops_tgfm_enrichment_summary.txt"), header=TRUE)
fig_4c <- mean_se_barplot_of_pops_score_binned_by_tgfm_pip(pops_summary_df, independent_traits) +theme(plot.margin = margin(5.5, 5.5, 0.0, 5.5, "points")) +
          theme(axis.text.x = element_text(angle = 35, hjust=1)) +
          labs(x = NULL)



# Fig 4D
ldl_enrichment_file <- paste0(ldl_silver_standard_gene_set_enrichment_dir, "ldl_silver_standard_enrichment_summary.txt")
fig_4d <- make_ldl_silver_standard_stratefied_gene_pip_emperical_fdr_plot(ldl_enrichment_file) +
		  labs(y="FDR") +
	      theme(plot.margin = margin(5.5, 5.5, 0.0, 5.5, "points")) 
fig_4d <- plot_grid(fig_4d, NULL, ncol=1,rel_heights=c(1,.247))


# MAke joint plot
fig_4ab = plot_grid(fig_4a, NULL, fig_4b, ncol=3, rel_widths=c(.6,.02,.45), labels=c("a","", "b"))
fig_4cd = plot_grid(fig_4c, NULL, fig_4d, ncol=3, rel_widths=c(.6,.02,.45), labels=c("c","", "d"))
fig_4 <- plot_grid(fig_4ab, fig_4cd, ncol=1, rel_heights=c(.6,.43))
output_file <- paste0(visualize_tgfm_dir, "figure4.pdf")
ggsave(fig_4, file=output_file, width=7.2, height=6.2, units="in")























