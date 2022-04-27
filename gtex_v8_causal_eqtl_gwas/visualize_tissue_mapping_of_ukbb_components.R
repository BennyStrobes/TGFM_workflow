args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(reshape2)
library(pheatmap)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

extract_component_pip_df <- function(studies, ukbb_susie_component_dir) {
	pip_arr <- c()
	study_arr <- c()
	for (study_iter in 1:length(studies)) {
		study = studies[study_iter]
		study_file <- paste0(ukbb_susie_component_dir, study, "_organized_susie_components.txt")
		study_df <- read.table(study_file, header=TRUE, sep="\t")

		pip_arr <- c(pip_arr, study_df$lead_variant_pi)
		study_arr <- c(study_arr, rep(study, length(study_df$lead_variant_pi)))
	}


	df <- data.frame(pip=pip_arr, study=factor(study_arr,levels=studies))
	return(df)	
}

extract_num_components_df <- function(studies, ukbb_susie_component_dir) {
	thresholds <- c(0.0, .1, .3, .5, .7, .9)
	num_components_arr <- c()
	study_arr <- c()
	threshold_arr <- c()
	for (study_iter in 1:length(studies)) {
		study = studies[study_iter]
		study_file <- paste0(ukbb_susie_component_dir, study, "_organized_susie_components.txt")
		study_df <- read.table(study_file, header=TRUE, sep="\t")
		for (thresh_iter in 1:length(thresholds)) {
			threshold <- thresholds[thresh_iter]
			num_components <- sum(study_df$lead_variant_pi > threshold)

			num_components_arr <- c(num_components_arr, num_components)
			study_arr <- c(study_arr, study)
			threshold_arr <- c(threshold_arr, threshold)
		}
	}


	df <- data.frame(num_components=num_components_arr, study=factor(study_arr,levels=studies), threshold=factor(threshold_arr, thresholds))
	return(df)
}

remove_roadmap_and_encode_labels_v2 <- function(df) {
	original_cts <- as.character(df$cell_type)
	print(original_cts)

	for (itera in 1:length(original_cts)) {
		ct <- original_cts[itera]
		if (endsWith(ct, ".Roadmap")) {
			new_ct <- substr(ct,1,nchar(ct)-8)
		}
		if (endsWith(ct, ".ENCODE")) {
			new_ct <- substr(ct,1,nchar(ct)-7)
		}
		if (new_ct == "A549_treated_with_ethanol_0.02_percent_for_1_hour") {
			new_ct = "A549"
		}
		if (new_ct == "H1_BMP4_Derived_Mesendoderm_Cultured_Cells") {
			new_ct = "H1_BMP4_mesendoderm"
		}
		if (new_ct == "H1_BMP4_Derived_Trophoblast_Cultured_Cells") {
			new_ct = "H1_BMP4_trophoblast"
		}
		if (new_ct == "H1_Derived_Mesenchymal_Stem_Cells") {
			new_ct = "H1_Mesenchymal"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct == "H1_Neuronal_Progenitor"
		}
		if (new_ct == "myotube_originated_from_skeletal_muscle_myoblast") {
			new_ct = "myotube"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct = "H1_Neuronal_Progenitor"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}
		if (new_ct == "CD8−positive_alpha−beta_T_cell") {
			new_ct = "CD8-positive_ab_T_cell"
		}
		if (new_ct == "endothelial_cell_of_umbilical_vein") {
			new_ct = "endothelial_umbilical_vein"
		}
		if (new_ct == "induced_pluripotent_stem_cell") {
			new_ct = "iPSC"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}

		original_cts[itera] = new_ct
	}
	df$cell_type = original_cts
	df$cell_type = factor(df$cell_type)

	return(df)
}


remove_roadmap_and_encode_labels <- function(df) {
	original_cts <- as.character(df$cell_type)
	print(original_cts)

	for (itera in 1:length(original_cts)) {
		ct <- original_cts[itera]
		if (endsWith(ct, "-Roadmap")) {
			new_ct <- substr(ct,1,nchar(ct)-8)
		}
		if (endsWith(ct, "-ENCODE")) {
			new_ct <- substr(ct,1,nchar(ct)-7)
		}
		if (new_ct == "A549_treated_with_ethanol_0.02_percent_for_1_hour") {
			new_ct = "A549"
		}
		if (new_ct == "H1_BMP4_Derived_Mesendoderm_Cultured_Cells") {
			new_ct = "H1_BMP4_mesendoderm"
		}
		if (new_ct == "H1_BMP4_Derived_Trophoblast_Cultured_Cells") {
			new_ct = "H1_BMP4_trophoblast"
		}
		if (new_ct == "H1_Derived_Mesenchymal_Stem_Cells") {
			new_ct = "H1_Mesenchymal"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct == "H1_Neuronal_Progenitor"
		}
		if (new_ct == "myotube_originated_from_skeletal_muscle_myoblast") {
			new_ct = "myotube"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct = "H1_Neuronal_Progenitor"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}
		if (new_ct == "CD8−positive_alpha−beta_T_cell") {
			new_ct = "CD8-positive_ab_T_cell"
		}
		if (new_ct == "endothelial_cell_of_umbilical_vein") {
			new_ct = "endothelial_umbilical_vein"
		}
		if (new_ct == "induced_pluripotent_stem_cell") {
			new_ct = "iPSC"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}

		original_cts[itera] = new_ct
	}
	df$cell_type = original_cts
	df$cell_type = factor(df$cell_type)

	return(df)
}

make_number_of_abc_links_per_cell_type_bar_plot <- function(num_epimap_links_df) {
	num_epimap_links_df <- num_epimap_links_df[num_epimap_links_df$threshold <= .1,]


	num_epimap_links_df <- remove_roadmap_and_encode_labels(num_epimap_links_df)

	p<-ggplot(data=num_epimap_links_df, aes(x=cell_type, y=num_links, fill=factor(threshold))) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="top") +
  			labs(y="Number of ABC links", fill="ABC threshold", x="ABC cell type") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8.5))	
  	return(p)	
}


make_number_of_epimap_links_per_cell_type_bar_plot <- function(num_epimap_links_df) {
	num_epimap_links_df <- num_epimap_links_df[num_epimap_links_df$threshold > .4,]

	p<-ggplot(data=num_epimap_links_df, aes(x=cell_type, y=num_links, fill=factor(threshold))) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="top") +
  			labs(y="Number of epimap links", fill="correlation threshold", x="Epimap cell type") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
  	return(p)	
}

make_number_of_components_per_study_bar_plot <- function(num_components_df) {
	p<-ggplot(data=num_components_df, aes(x=study, y=num_components, fill=threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="top") +
  			labs(y="Number of SuSiE components", fill="PIP threshold", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
  	return(p)
}

make_pip_per_susie_component_boxplot <- function(component_pip_df) {
	p<-ggplot(data=component_pip_df, aes(x=study, y=pip)) +
  			geom_boxplot() +
  			figure_theme() + 
  			theme(legend.position="top") +
  			labs(y="SuSiE PIP", x="study") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
  	return(p)	
}

extract_num_links_df <- function(file_name) {
	df <- read.table(file_name, header=TRUE, sep="\t")
	return(df)
}


make_fraction_of_components_with_cell_type_support_bar_plot_in_single_study_vary_thresh <- function(study_file,study) {
	thresholds <- c(0, .2, .4)


	total_weights_arr <- c()
	ct_arr <- c()
	threshold_arr <- c()

	for (threshold_iter in 1:length(thresholds)) {
	threshold <- thresholds[threshold_iter]

	df <- read.table(study_file, header=TRUE,sep="\t")

	df <- df[df$max_un_normalized_cell_type_weight >= threshold, ]
	num_col <- dim(df)[2]
	num_components <- dim(df)[1]
	df <- df[,9:num_col]
	melted_df <- melt(df)
	cell_types <- as.character(sort(unique(melted_df$variable)))

	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		weights <- melted_df$value[as.character(melted_df$variable)==cell_type]
		total_weights_arr <- c(total_weights_arr, sum(weights)/num_components)
		ct_arr <- c(ct_arr, cell_type)
		threshold_arr <- c(threshold_arr, threshold)
	}
	}
	final_df <- data.frame(fraction_causal=total_weights_arr, cell_type=factor(ct_arr), threshold=factor(threshold_arr, thresholds))

	p<-ggplot(data=final_df, aes(x=cell_type, y=fraction_causal, fill=threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Fraction of\nSuSiE components causal", x="Epimap cell type", title=study) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
  	return(p)	
}



make_fraction_of_components_with_cell_type_support_bar_plot_in_single_study <- function(study_file,study) {
	df <- read.table(study_file, header=TRUE,sep="\t")
	num_col <- dim(df)[2]
	num_components <- dim(df)[1]
	df <- df[,9:num_col]
	melted_df <- melt(df)
	cell_types <- as.character(sort(unique(melted_df$variable)))

	total_weights_arr <- c()
	ct_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		weights <- melted_df$value[as.character(melted_df$variable)==cell_type]
		total_weights_arr <- c(total_weights_arr, sum(weights)/num_components)
		ct_arr <- c(ct_arr, cell_type)
	}
	final_df <- data.frame(fraction_causal=total_weights_arr, cell_type=factor(ct_arr))

	p<-ggplot(data=final_df, aes(x=cell_type, y=fraction_causal)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Fraction of\nSuSiE components causal", x="Epimap cell type", title=study) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
  	return(p)	
}

make_fraction_of_components_with_abc_cell_type_support_bar_plot_in_single_study_v2 <- function(study_file,study) {
	df <- read.table(study_file, header=TRUE,sep="\t")
	num_col <- dim(df)[2]
	num_components <- dim(df)[1]
	df <- df[,9:num_col]
	melted_df <- melt(df)
	cell_types <- as.character(sort(unique(melted_df$variable)))

	total_weights_arr <- c()
	ct_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		weights <- melted_df$value[as.character(melted_df$variable)==cell_type]
		total_weights_arr <- c(total_weights_arr, sum(weights)/num_components)
		ct_arr <- c(ct_arr, cell_type)
	}
	final_df <- data.frame(fraction_causal=total_weights_arr, cell_type=factor(ct_arr, levels=ct_arr[order(total_weights_arr)]))


	final_df <- remove_roadmap_and_encode_labels_v2(final_df)

	ct_final = as.character(final_df$cell_type)

	final_df2 <- data.frame(fraction_causal=final_df$fraction_causal, cell_type=factor(ct_final, levels=ct_final[order(-final_df$fraction_causal)]))

	p<-ggplot(data=final_df2, aes(x=cell_type, y=fraction_causal)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Fraction of\nSuSiE components causal", x="ABC cell type", title=study) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
  	return(p)	
}

make_fraction_of_components_with_abc_cell_type_support_bar_plot_in_single_study <- function(study_file,study) {
	df <- read.table(study_file, header=TRUE,sep="\t")
	num_col <- dim(df)[2]
	num_components <- dim(df)[1]
	df <- df[,9:num_col]
	melted_df <- melt(df)
	cell_types <- as.character(sort(unique(melted_df$variable)))

	total_weights_arr <- c()
	ct_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		weights <- melted_df$value[as.character(melted_df$variable)==cell_type]
		total_weights_arr <- c(total_weights_arr, sum(weights)/num_components)
		ct_arr <- c(ct_arr, cell_type)
	}
	final_df <- data.frame(fraction_causal=total_weights_arr, cell_type=factor(ct_arr))

	final_df <- remove_roadmap_and_encode_labels_v2(final_df)

	p<-ggplot(data=final_df, aes(x=cell_type, y=fraction_causal)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Fraction of\nSuSiE components causal", x="ABC cell type", title=study) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))	
  	return(p)	
}


make_fraction_of_components_with_cell_type_support_vs_num_epimap_components_scatter_in_single_study <- function(study_file,study, num_epimap_links_df) {
	num_epimap_links_df <- num_epimap_links_df[num_epimap_links_df$threshold==.5,]
	df <- read.table(study_file, header=TRUE,sep="\t")
	num_col <- dim(df)[2]
	num_components <- dim(df)[1]
	df <- df[,9:num_col]
	melted_df <- melt(df)
	cell_types <- as.character(sort(unique(melted_df$variable)))

	total_weights_arr <- c()
	num_components_arr <- c()
	ct_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		weights <- melted_df$value[as.character(melted_df$variable)==cell_type]
		total_weights_arr <- c(total_weights_arr, sum(weights)/num_components)
		ct_arr <- c(ct_arr, cell_type)

		value_arr = num_epimap_links_df$num_links[num_epimap_links_df$cell_type==cell_type]
		num_components_arr <- c(num_components_arr, value_arr[1])
	}
	final_df <- data.frame(fraction_causal=total_weights_arr, num_components=num_components_arr, cell_type=factor(ct_arr))	


	p <- ggplot(final_df, aes(x=num_components, y=fraction_causal)) +
 		 geom_point() +
 		 figure_theme() +
 		 labs(x="Number of epimap links (in a cell type)", y="Fraction components epimap causal\n(in a cell type)", title=study)

 	return(p)
}




make_fraction_of_components_with_cell_type_support_vs_num_epimap_components_scatter <- function(studies, epimap_ukbb_overlap_dir, num_epimap_links_df) {

	plot_list <- list()

	for (study_iter in 1:length(studies)) {
		study <- studies[study_iter]
		study_file <- paste0(epimap_ukbb_overlap_dir, study, "_component_overlap_with_epimap_0.75.txt")
		p <- make_fraction_of_components_with_cell_type_support_vs_num_epimap_components_scatter_in_single_study(study_file,study, num_epimap_links_df)
		plot_list[[study_iter]] = p
	}

	merged_plot <- plot_grid(plotlist=plot_list, ncol=1)

	return(merged_plot)	
}


make_fraction_of_components_with_cell_type_support_bar_plot_vary_thresh <- function(studies, epimap_ukbb_overlap_dir) {
	# Temper

	plot_list <- list()

	for (study_iter in 1:length(studies)) {
		study <- studies[study_iter]
		study_file <- paste0(epimap_ukbb_overlap_dir, study, "_component_overlap_with_epimap_0.75.txt")
		p <- make_fraction_of_components_with_cell_type_support_bar_plot_in_single_study_vary_thresh(study_file,study)
		plot_list[[study_iter]] = p
	}

	merged_plot <- plot_grid(plotlist=plot_list, ncol=1)

	return(merged_plot)


}

make_fraction_of_components_with_cell_type_support_bar_plot <- function(studies, epimap_ukbb_overlap_dir) {
	# Temper

	plot_list <- list()

	for (study_iter in 1:length(studies)) {
		study <- studies[study_iter]
		study_file <- paste0(epimap_ukbb_overlap_dir, study, "_component_overlap_with_epimap_0.75.txt")
		p <- make_fraction_of_components_with_cell_type_support_bar_plot_in_single_study(study_file,study)
		plot_list[[study_iter]] = p
	}

	merged_plot <- plot_grid(plotlist=plot_list, ncol=1)

	return(merged_plot)


}

make_fraction_of_components_with_abc_cell_type_support_bar_plot <- function(studies, epimap_ukbb_overlap_dir) {
	# Temper

	plot_list <- list()

	for (study_iter in 1:length(studies)) {
		study <- studies[study_iter]
		study_file <- paste0(epimap_ukbb_overlap_dir, study, "_component_overlap_with_abc_0.015.txt")
		p <- make_fraction_of_components_with_abc_cell_type_support_bar_plot_in_single_study_v2(study_file,study)
		plot_list[[study_iter]] = p
	}

	merged_plot <- plot_grid(plotlist=plot_list, ncol=1)

	return(merged_plot)


}

compare_fine_mapped_vs_sldsc_abc_enrichments <- function(abc_paper_study_name, abc_paper_enrichment_file) {
	print(abc_paper_enrichment_file)
	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)

	abc_paper_df_subset <- abc_paper_df[as.character(abc_paper_df$trait) == abc_paper_study_name,]

	print(cor(abc_paper_df_subset$FM.Enrichment, abc_paper_df_subset$LDSC.Enrichment))
	print(cor(abc_paper_df_subset$FM.log10p, log10(abc_paper_df_subset$LDSC.EnrichmentP)))
}

compare_raw_abc_enrichments <- function(my_study_name, abc_paper_study_name, abc_ukbb_overlap_dir, abc_paper_enrichment_file) {
	my_abc_file <- paste0(abc_ukbb_overlap_dir, my_study_name, "_raw_component_overlap_with_abc_0.015.txt")
	
	df <- read.table(my_abc_file, header=TRUE,sep="\t")
	num_col <- dim(df)[2]
	num_components <- dim(df)[1]
	pis <- as.character(df$cs_pi)
	num_variants_per_component <- c()

	for (component_num in 1:num_components) {
		component_pi <- pis[component_num]
		num_variants <- length(strsplit(component_pi, ",")[[1]])
		num_variants_per_component <- c(num_variants_per_component, num_variants)
	}


	df <- df[,8:num_col]
	num_col <- dim(df)[2]
	for (col_num in 1:num_col) {
		df[, col_num] = df[, col_num]
	}

	melted_df <- melt(df)
	cell_types <- as.character(sort(unique(melted_df$variable)))


	total_weights_arr <- c()
	ct_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		weights <- melted_df$value[as.character(melted_df$variable)==cell_type]
		total_weights_arr <- c(total_weights_arr, sum(weights)/sum(num_variants_per_component))
		ct_arr <- c(ct_arr, cell_type)
	}
	final_df <- data.frame(fraction_causal=total_weights_arr, cell_type=factor(ct_arr))


	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)

	abc_paper_df_subset <- abc_paper_df[as.character(abc_paper_df$trait) == abc_paper_study_name,]

	abc_paper_df_subset$Biosample = str_replace_all(as.character(abc_paper_df_subset$Biosample), "-", ".")

	my_ratio_arr <- c()
	abc_paper_ratio_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		my_ratio <- total_weights_arr[cell_type_iter]
		cell_type <- cell_types[cell_type_iter]

		abc_paper_indices <- abc_paper_df_subset$Biosample == cell_type

		if (sum(abc_paper_indices) == 1.0) {
			abc_paper_overlaps = abc_paper_df_subset$FM.nVariantOverlap[abc_paper_indices]
			abc_paper_variants = abc_paper_df_subset$FM.nVariants[abc_paper_indices]
			abc_paper_ratio <- abc_paper_overlaps/abc_paper_variants
			abc_paper_ratio <- abc_paper_df_subset$FM.Enrichment[abc_paper_indices]
			my_ratio_arr <- c(my_ratio_arr, my_ratio)
			abc_paper_ratio_arr <- c(abc_paper_ratio_arr, abc_paper_ratio)
		}
	}
	print(cor(my_ratio_arr, abc_paper_ratio_arr))
}



compare_abc_enrichments <- function(my_study_name, abc_paper_study_name, abc_ukbb_overlap_dir, abc_paper_enrichment_file) {
	my_abc_file <- paste0(abc_ukbb_overlap_dir, my_study_name, "_component_overlap_with_abc_0.015.txt")
	df <- read.table(my_abc_file, header=TRUE,sep="\t")
	df <- df[df$lead_variant_pi > .1,]
	num_col <- dim(df)[2]
	num_components <- dim(df)[1]
	weights = df$sum_un_normalized_cell_type_weight

	df <- df[,9:num_col]
	num_col <- dim(df)[2]
	for (col_num in 1:num_col) {
		df[, col_num] = df[, col_num]*weights
	}

	melted_df <- melt(df)
	cell_types <- as.character(sort(unique(melted_df$variable)))


	total_weights_arr <- c()
	ct_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		weights <- melted_df$value[as.character(melted_df$variable)==cell_type]
		total_weights_arr <- c(total_weights_arr, sum(weights)/num_components)
		ct_arr <- c(ct_arr, cell_type)
	}
	final_df <- data.frame(fraction_causal=total_weights_arr, cell_type=factor(ct_arr))


	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)


	#abc_paper_df_subset <- abc_paper_df[as.character(abc_paper_df$trait) == "IBD",]
	#print(abc_paper_df_subset)


	abc_paper_df_subset <- abc_paper_df[as.character(abc_paper_df$trait) == abc_paper_study_name,]

	abc_paper_df_subset$Biosample = str_replace_all(as.character(abc_paper_df_subset$Biosample), "-", ".")

	my_ratio_arr <- c()
	abc_paper_ratio_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		my_ratio <- total_weights_arr[cell_type_iter]
		cell_type <- cell_types[cell_type_iter]

		abc_paper_indices <- abc_paper_df_subset$Biosample == cell_type

		if (sum(abc_paper_indices) == 1.0) {
			abc_paper_overlaps = abc_paper_df_subset$FM.nVariantOverlap[abc_paper_indices]
			abc_paper_variants = abc_paper_df_subset$FM.nVariants[abc_paper_indices]
			abc_paper_ratio <- abc_paper_overlaps/abc_paper_variants
			#abc_paper_ratio <- abc_paper_df_subset$FM.Enrichment[abc_paper_indices]
			my_ratio_arr <- c(my_ratio_arr, my_ratio)
			abc_paper_ratio_arr <- c(abc_paper_ratio_arr, abc_paper_ratio)
		}
	}

	corry <- cor(my_ratio_arr, abc_paper_ratio_arr)
	plot_df <- data.frame(my_ratio=my_ratio_arr, abc_paper_ratio=abc_paper_ratio_arr)

	p <- ggplot(plot_df, aes(x=my_ratio, y=abc_paper_ratio)) +
 		 geom_point() +
 		 figure_theme() +
 		 geom_abline() +
 		 labs(x="Fraction causal", y="Fraction causal (ABC paper)", title=my_study_name)

 	return(p)

}

paper_sldsc_enrichment_heatmap <- function(abc_paper_enrichment_file, output_file_name) {
	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)


	final_df <- data.frame(trait=as.character(abc_paper_df$trait), biosample=as.character(abc_paper_df$Biosample), enrichment=abc_paper_df$LDSC.Enrichment)
	#print(matrix(final_df))

	original_mat <- dcast(final_df, trait ~ biosample)

	matty <- original_mat[, 2:dim(original_mat)[2]]
	rownames(matty) = original_mat$trait
	valid_cols <- c()
	column_names <- colnames(matty)
	for (itera in 1:length(column_names)) {
		column_name <- column_names[itera]
		if (endsWith(column_name, "Roadmap")) {
			valid_cols <- c(valid_cols, TRUE)
		} else if (endsWith(column_name, "ENCODE")) {
			valid_cols <- c(valid_cols, TRUE)
		} else {
			valid_cols <- c(valid_cols, FALSE)
		}
	}

	row_names <- rownames(matty)
	valid_rows <- row_names != "CRC"
	matty <- matty[valid_rows,]

	matty <- matty[, valid_cols]
	column_names <- colnames(matty)


	for (itera in 1:length(column_names)) {
		ct <- column_names[itera]
		if (endsWith(ct, "Roadmap")) {
			new_ct <- substr(ct,1,nchar(ct)-8)
		}
		else if (endsWith(ct, "ENCODE")) {
			new_ct <- substr(ct,1,nchar(ct)-7)
		} else {
			new_ct = ct
		}
		column_names[itera] = new_ct
	}
	colnames(matty) <- column_names

	#print(head(as.numeric(matty)))
	#print(head(as.matrix(matty)))
	pdf(file=output_file_name, width=12, height=13)
	p <- pheatmap((as.matrix(matty)), scale = "none")
	dev.off()
}
paper_fine_mapped_fraction_linked_heatmap <- function(abc_paper_enrichment_file, output_file_name) {
	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)

	final_df <- data.frame(trait=as.character(abc_paper_df$trait), biosample=as.character(abc_paper_df$Biosample), enrichment=abc_paper_df$FM.nVariantOverlaps/abc_paper_df$FM.nVariants)
	#print(matrix(final_df))

	original_mat <- dcast(final_df, trait ~ biosample)

	matty <- original_mat[, 2:dim(original_mat)[2]]
	rownames(matty) = original_mat$trait
	valid_cols <- c()
	column_names <- colnames(matty)
	for (itera in 1:length(column_names)) {
		column_name <- column_names[itera]
		if (endsWith(column_name, "Roadmap")) {
			valid_cols <- c(valid_cols, TRUE)
		} else if (endsWith(column_name, "ENCODE")) {
			valid_cols <- c(valid_cols, TRUE)
		} else {
			valid_cols <- c(valid_cols, FALSE)
		}
	}

	row_names <- rownames(matty)
	valid_rows <- row_names != "CRC"
	matty <- matty[valid_rows,]

	matty <- matty[, valid_cols]
	column_names <- colnames(matty)


	for (itera in 1:length(column_names)) {
		ct <- column_names[itera]
		if (endsWith(ct, "Roadmap")) {
			new_ct <- substr(ct,1,nchar(ct)-8)
		}
		else if (endsWith(ct, "ENCODE")) {
			new_ct <- substr(ct,1,nchar(ct)-7)
		} else {
			new_ct = ct
		}
		column_names[itera] = new_ct
	}
	colnames(matty) <- column_names

	#print(head(as.numeric(matty)))
	#print(head(as.matrix(matty)))
	pdf(file=output_file_name, width=12, height=13)
	p <- pheatmap((as.matrix(matty)), scale = "none")
	dev.off()
}

abc_non_causal_overlap_barplot <- function(abc_paper_study_name, my_study_name, abc_paper_enrichment_file) {
	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)


	#df <- data.frame(trait=as.character(abc_paper_df$trait), biosample=as.character(abc_paper_df$Biosample), enrichment=abc_paper_df$FM.Enrichment)	
	df <- data.frame(trait=as.character(abc_paper_df$trait), biosample=as.character(abc_paper_df$Biosample), enrichment=abc_paper_df$FM.nVariantOverlap/abc_paper_df$FM.nVariants)	




	df <- df[as.character(df$trait) == abc_paper_study_name,]

	biosamples <- as.character(df$biosample)
	valid_rows <- c()

	for (row_num in 1:length(biosamples)) {
		biosample_name <- biosamples[row_num]
		if (endsWith(biosample_name, "Roadmap")) {
			valid_rows <- c(valid_rows, TRUE)
		} else if (endsWith(biosample_name, "ENCODE")) {
			valid_rows <- c(valid_rows, TRUE)
		} else {
			valid_rows <- c(valid_rows, FALSE)
		}		
	}

	df <- df[valid_rows,]

	biosamples <- as.character(df$biosample)
	new_biosample_names <- c()

	for (itera in 1:length(biosamples)) {
		ct <- biosamples[itera]
		if (endsWith(ct, "Roadmap")) {
			new_ct <- substr(ct,1,nchar(ct)-8)
		}
		else if (endsWith(ct, "ENCODE")) {
			new_ct <- substr(ct,1,nchar(ct)-7)
		} else {
			new_ct = ct
		}


		if (new_ct == "A549_treated_with_ethanol_0.02_percent_for_1_hour") {
			new_ct = "A549"
		}
		if (new_ct == "H1_BMP4_Derived_Mesendoderm_Cultured_Cells") {
			new_ct = "H1_BMP4_mesendoderm"
		}
		if (new_ct == "H1_BMP4_Derived_Trophoblast_Cultured_Cells") {
			new_ct = "H1_BMP4_trophoblast"
		}
		if (new_ct == "H1_Derived_Mesenchymal_Stem_Cells") {
			new_ct = "H1_Mesenchymal"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct == "H1_Neuronal_Progenitor"
		}
		if (new_ct == "myotube_originated_from_skeletal_muscle_myoblast") {
			new_ct = "myotube"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct = "H1_Neuronal_Progenitor"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}
		if (new_ct == "CD8−positive_alpha−beta_T_cell") {
			new_ct = "CD8-positive_ab_T_cell"
		}
		if (new_ct == "endothelial_cell_of_umbilical_vein") {
			new_ct = "endothelial_umbilical_vein"
		}
		if (new_ct == "induced_pluripotent_stem_cell") {
			new_ct = "iPSC"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}



		new_biosample_names[itera] = new_ct
	}
	
	df$biosample <- as.character(new_biosample_names)

	final_df <- data.frame(enrichment=df$enrichment, cell_type=factor(df$biosample, levels=df$biosample[order(-df$enrichment)]))

	p<-ggplot(data=final_df, aes(x=cell_type, y=enrichment)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Non Causal overlap", x="ABC cell type", title=my_study_name) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8))
  	return(p)	

}




abc_enrichment_barplot <- function(abc_paper_study_name, my_study_name, abc_paper_enrichment_file) {
	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)


	df <- data.frame(trait=as.character(abc_paper_df$trait), biosample=as.character(abc_paper_df$Biosample), enrichment=abc_paper_df$FM.Enrichment)	
	#df <- data.frame(trait=as.character(abc_paper_df$trait), biosample=as.character(abc_paper_df$Biosample), enrichment=abc_paper_df$FM.nVariantOverlap/abc_paper_df$FM.nVariants)	




	df <- df[as.character(df$trait) == abc_paper_study_name,]

	biosamples <- as.character(df$biosample)
	valid_rows <- c()

	for (row_num in 1:length(biosamples)) {
		biosample_name <- biosamples[row_num]
		if (endsWith(biosample_name, "Roadmap")) {
			valid_rows <- c(valid_rows, TRUE)
		} else if (endsWith(biosample_name, "ENCODE")) {
			valid_rows <- c(valid_rows, TRUE)
		} else {
			valid_rows <- c(valid_rows, FALSE)
		}		
	}

	df <- df[valid_rows,]

	biosamples <- as.character(df$biosample)
	new_biosample_names <- c()

	for (itera in 1:length(biosamples)) {
		ct <- biosamples[itera]
		if (endsWith(ct, "Roadmap")) {
			new_ct <- substr(ct,1,nchar(ct)-8)
		}
		else if (endsWith(ct, "ENCODE")) {
			new_ct <- substr(ct,1,nchar(ct)-7)
		} else {
			new_ct = ct
		}


		if (new_ct == "A549_treated_with_ethanol_0.02_percent_for_1_hour") {
			new_ct = "A549"
		}
		if (new_ct == "H1_BMP4_Derived_Mesendoderm_Cultured_Cells") {
			new_ct = "H1_BMP4_mesendoderm"
		}
		if (new_ct == "H1_BMP4_Derived_Trophoblast_Cultured_Cells") {
			new_ct = "H1_BMP4_trophoblast"
		}
		if (new_ct == "H1_Derived_Mesenchymal_Stem_Cells") {
			new_ct = "H1_Mesenchymal"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct == "H1_Neuronal_Progenitor"
		}
		if (new_ct == "myotube_originated_from_skeletal_muscle_myoblast") {
			new_ct = "myotube"
		}
		if (new_ct == "H1_Derived_Neuronal_Progenitor_Cultured_Cells") {
			new_ct = "H1_Neuronal_Progenitor"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}
		if (new_ct == "CD8−positive_alpha−beta_T_cell") {
			new_ct = "CD8-positive_ab_T_cell"
		}
		if (new_ct == "endothelial_cell_of_umbilical_vein") {
			new_ct = "endothelial_umbilical_vein"
		}
		if (new_ct == "induced_pluripotent_stem_cell") {
			new_ct = "iPSC"
		}
		if (new_ct == "CD56−positive_natural_killer_cells") {
			new_ct = "CD56-positive_NK"
		}



		new_biosample_names[itera] = new_ct
	}
	
	df$biosample <- as.character(new_biosample_names)

	final_df <- data.frame(enrichment=df$enrichment, cell_type=factor(df$biosample, levels=df$biosample[order(-df$enrichment)]))

	p<-ggplot(data=final_df, aes(x=cell_type, y=enrichment)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Enrichment", x="ABC cell type", title=my_study_name) +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8))	+
  			geom_hline(yintercept=1.0,color='red')
  	return(p)	

}


paper_fine_mapped_enrichment_heatmap <- function(abc_paper_enrichment_file, output_file_name) {
	abc_paper_df <- read.table(abc_paper_enrichment_file, header=TRUE)

	final_df <- data.frame(trait=as.character(abc_paper_df$trait), biosample=as.character(abc_paper_df$Biosample), enrichment=abc_paper_df$FM.Enrichment)
	#print(matrix(final_df))

	original_mat <- dcast(final_df, trait ~ biosample)

	matty <- original_mat[, 2:dim(original_mat)[2]]
	rownames(matty) = original_mat$trait
	valid_cols <- c()
	column_names <- colnames(matty)
	for (itera in 1:length(column_names)) {
		column_name <- column_names[itera]
		if (endsWith(column_name, "Roadmap")) {
			valid_cols <- c(valid_cols, TRUE)
		} else if (endsWith(column_name, "ENCODE")) {
			valid_cols <- c(valid_cols, TRUE)
		} else {
			valid_cols <- c(valid_cols, FALSE)
		}
	}

	row_names <- rownames(matty)
	valid_rows <- row_names != "CRC"
	matty <- matty[valid_rows,]

	matty <- matty[, valid_cols]
	column_names <- colnames(matty)


	for (itera in 1:length(column_names)) {
		ct <- column_names[itera]
		if (endsWith(ct, "Roadmap")) {
			new_ct <- substr(ct,1,nchar(ct)-8)
		}
		else if (endsWith(ct, "ENCODE")) {
			new_ct <- substr(ct,1,nchar(ct)-7)
		} else {
			new_ct = ct
		}
		column_names[itera] = new_ct
	}
	colnames(matty) <- column_names

	#print(head(as.numeric(matty)))
	#print(head(as.matrix(matty)))
	pdf(file=output_file_name, width=12, height=13)
	p <- pheatmap((as.matrix(matty)), scale = "none")
	dev.off()
}


output_dir <- args[1]
ukbb_susie_component_dir <- args[2]
epimap_ukbb_overlap_dir <- args[3]
abc_ukbb_overlap_dir <- args[4]
abc_paper_enrichment_file <- args[5]



# UKBB GWAS Studies
studies <- c("blood_WHITE_COUNT", "body_WHRadjBMIz", "bp_DIASTOLICadjMEDz", "lung_FEV1FVCzSMOKE", "repro_ChildLess")


# Extract df containing number of components per study
num_components_df <- extract_num_components_df(studies, ukbb_susie_component_dir)

# Extract df containing pips per study
component_pip_df <- extract_component_pip_df(studies, ukbb_susie_component_dir)

# Extract df containing num epimap links
num_epimap_links_df <- extract_num_links_df(paste0(epimap_ukbb_overlap_dir, "number_epimap_links.txt"))

# Extract df containing num abc links
num_abc_links_df <- extract_num_links_df(paste0(abc_ukbb_overlap_dir, "number_abc_links.txt"))


if (FALSE) {

# Visualize number of components per study
output_file <- paste0(output_dir, "number_of_susie_components_per_study.pdf")
bar_plot <- make_number_of_components_per_study_bar_plot(num_components_df)
ggsave(bar_plot, file=output_file, width=7.2, height=5.0, units="in")


# Visualize PIP boxplot
output_file <- paste0(output_dir, "pip_per_susie_component_boxplot.pdf")
boxplot <- make_pip_per_susie_component_boxplot(component_pip_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")


# Make number of abc links bar plot
output_file <- paste0(output_dir, "number_of_abc_links_per_cell_type.pdf")
bar_plot <- make_number_of_abc_links_per_cell_type_bar_plot(num_abc_links_df)
ggsave(bar_plot, file=output_file, width=10.2, height=5.0, units="in")
}

# Fraction of components with cell type specific ABC support
output_file <- paste0(output_dir, "fraction_of_components_with_abc_cell_type_support.pdf")
bar_plot <- make_fraction_of_components_with_abc_cell_type_support_bar_plot(studies[1:(length(studies)-1)], abc_ukbb_overlap_dir)
ggsave(bar_plot, file=output_file, width=10.2, height=20.3, units="in")

# Compare my ABC enrichments with their abc enrichments
my_study_name <- "body_WHRadjBMIz"
abc_paper_study_name <- "WHRadjBMI"
output_file <- paste0(output_dir, "compare_fraction_causal_with_abc_paper_", my_study_name, ".pdf")
scatter <- compare_abc_enrichments(my_study_name, abc_paper_study_name, abc_ukbb_overlap_dir, abc_paper_enrichment_file)
ggsave(scatter, file=output_file, width=7.2, height=5.0, units="in")

my_study_name <- "bp_DIASTOLICadjMEDz"
abc_paper_study_name <- "DBP"
output_file <- paste0(output_dir, "compare_fraction_causal_with_abc_paper_", my_study_name, ".pdf")
scatter <- compare_abc_enrichments(my_study_name, abc_paper_study_name, abc_ukbb_overlap_dir, abc_paper_enrichment_file)
ggsave(scatter, file=output_file, width=7.2, height=5.0, units="in")

my_study_name <- "blood_WHITE_COUNT"
abc_paper_study_name <- "WBC"
output_file <- paste0(output_dir, "compare_fraction_causal_with_abc_paper_", my_study_name, ".pdf")
scatter <- compare_abc_enrichments(my_study_name, abc_paper_study_name, abc_ukbb_overlap_dir, abc_paper_enrichment_file)
ggsave(scatter, file=output_file, width=7.2, height=5.0, units="in")

my_study_name <- "lung_FEV1FVCzSMOKE"
abc_paper_study_name <- "FEV1FVC"
output_file <- paste0(output_dir, "compare_fraction_causal_with_abc_paper_", my_study_name, ".pdf")
scatter <- compare_abc_enrichments(my_study_name, abc_paper_study_name, abc_ukbb_overlap_dir, abc_paper_enrichment_file)
ggsave(scatter, file=output_file, width=7.2, height=5.0, units="in")




# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "WHRadjBMI"
my_study_name <- "body_WHRadjBMIz"
output_file <- paste0(output_dir, "enrichment_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_enrichment_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")

# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "WBC"
my_study_name <- "blood_WHITE_COUNT"
output_file <- paste0(output_dir, "enrichment_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_enrichment_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")

# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "DBP"
my_study_name <- "bp_DIASTOLICadjMEDz"
output_file <- paste0(output_dir, "enrichment_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_enrichment_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")

# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "FEV1FVC"
my_study_name <- "lung_FEV1FVCzSMOKE"
output_file <- paste0(output_dir, "enrichment_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_enrichment_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")



# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "WHRadjBMI"
my_study_name <- "body_WHRadjBMIz"
output_file <- paste0(output_dir, "non_causal_overlap_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_non_causal_overlap_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")

# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "WBC"
my_study_name <- "blood_WHITE_COUNT"
output_file <- paste0(output_dir, "non_causal_overlap_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_non_causal_overlap_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")

# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "DBP"
my_study_name <- "bp_DIASTOLICadjMEDz"
output_file <- paste0(output_dir, "non_causal_overlap_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_non_causal_overlap_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")

# Compare my ABC enrichments with their abc enrichments
abc_paper_study_name <- "FEV1FVC"
my_study_name <- "lung_FEV1FVCzSMOKE"
output_file <- paste0(output_dir, "non_causal_overlap_of_components_with_abc_support_", my_study_name, ".pdf")
batplot <- abc_non_causal_overlap_barplot(abc_paper_study_name, my_study_name, abc_paper_enrichment_file)
ggsave(batplot, file=output_file, width=7.2, height=5.0, units="in")





# fine-mapped fraction linked heatmap
output_file <- paste0(output_dir, "abc_paper_fine_mapped_fraction_linked_heatmap.pdf")
paper_fine_mapped_fraction_linked_heatmap(abc_paper_enrichment_file, output_file)

# fine-mapped enrichment heatmap
output_file <- paste0(output_dir, "abc_paper_fine_mapped_enrichment_heatmap.pdf")
paper_fine_mapped_enrichment_heatmap(abc_paper_enrichment_file, output_file)

# fine-mapped enrichment heatmap
output_file <- paste0(output_dir, "abc_paper_sldsc_enrichment_heatmap.pdf")
paper_sldsc_enrichment_heatmap(abc_paper_enrichment_file, output_file)



if (FALSE) {

# Make number of epimap links bar plot
output_file <- paste0(output_dir, "number_of_epimap_links_per_cell_type.pdf")
bar_plot <- make_number_of_epimap_links_per_cell_type_bar_plot(num_epimap_links_df)
ggsave(bar_plot, file=output_file, width=7.2, height=5.0, units="in")


# Fraction of components with cell type specific epimap support
output_file <- paste0(output_dir, "fraction_of_components_with_epimap_cell_type_support.pdf")
bar_plot <- make_fraction_of_components_with_cell_type_support_bar_plot(studies[1:(length(studies)-1)], epimap_ukbb_overlap_dir)
ggsave(bar_plot, file=output_file, width=7.2, height=9.3, units="in")


# Fraction of components with cell type specific epimap support
output_file <- paste0(output_dir, "fraction_of_components_with_epimap_cell_type_support_vary_thresh.pdf")
bar_plot <- make_fraction_of_components_with_cell_type_support_bar_plot_vary_thresh(studies[1:(length(studies)-1)], epimap_ukbb_overlap_dir)
ggsave(bar_plot, file=output_file, width=7.2, height=9.3, units="in")


# Fraction of components with cell type specific epimap support vs number epimap components scatter
output_file <- paste0(output_dir, "fraction_of_components_with_epimap_cell_type_support_vs_number_epimap_components_scatter.pdf")
bar_plot <- make_fraction_of_components_with_cell_type_support_vs_num_epimap_components_scatter(studies[1:(length(studies)-1)], epimap_ukbb_overlap_dir, num_epimap_links_df)
ggsave(bar_plot, file=output_file, width=7.2, height=9.3, units="in")

}
