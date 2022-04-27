args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(hash)
library(PRROC)
library(pROC)
library(dplyr)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



make_boxplot_comparing_nominal_log10_pvalues_depending_on_whether_there_was_a_susie_component <- function(ce_df) {
	pvalue_vec <- c()
	susie_component_vec <- c()

	num_samples <- dim(ce_df)[1]

	indices <- (ce_df$num_trait_components == 0) & (ce_df$num_tissue_components == 0)
	temp_pvalue <- ce_df$nominal_pvalue[indices]
	temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]
	pvalue_vec <- c(pvalue_vec, temp_pvalue)
	susie_component_vec <- c(susie_component_vec, rep("no_component", length(temp_pvalue)))

	indices <- (ce_df$num_trait_components != 0) & (ce_df$num_tissue_components == 0)
	temp_pvalue <- ce_df$nominal_pvalue[indices]
	temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]
	pvalue_vec <- c(pvalue_vec, temp_pvalue)
	susie_component_vec <- c(susie_component_vec, rep("trait_component_only", length(temp_pvalue)))

	indices <- (ce_df$num_trait_components == 0) & (ce_df$num_tissue_components != 0)
	temp_pvalue <- ce_df$nominal_pvalue[indices]
	temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]
	pvalue_vec <- c(pvalue_vec, temp_pvalue)
	susie_component_vec <- c(susie_component_vec, rep("eqtl_component_only", length(temp_pvalue)))


	indices <- (ce_df$num_trait_components != 0) & (ce_df$num_tissue_components != 0)
	temp_pvalue <- ce_df$nominal_pvalue[indices]
	temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]
	pvalue_vec <- c(pvalue_vec, temp_pvalue)
	susie_component_vec <- c(susie_component_vec, rep("trait_eqtl_component", length(temp_pvalue)))


	df <- data.frame(pvalue=-log10(pvalue_vec+1e-30), component_status=factor(susie_component_vec, levels=c("no_component", "eqtl_component_only", "trait_component_only", "trait_eqtl_component")))

	p <- ggplot(data=df, aes(x=component_status, y=pvalue, fill=component_status)) +
		geom_boxplot() +
		figure_theme() +
		theme(legend.position="none") + 
		labs(y="-log10(pvalue) for (trait, gene, tissue)", fill="", x="") +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
	return(p)

}

make_boxplot_comparing_nominal_log10_pvalues_stratefied_by_tissue_type <- function(ce_df) {
	pvalue_vec <- c()
	tissue_vec <- c()

	num_samples <- dim(ce_df)[1]

	tissue_types <- as.character(sort(unique(ce_df$tissue)))


	for (tissue_iter in 1:length(tissue_types)) {
		tissue_type <- tissue_types[tissue_iter]

		indices <- ce_df$tissue == tissue_type
		temp_pvalue <- ce_df$nominal_pvalue[indices]
		temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]
		pvalue_vec <- c(pvalue_vec, temp_pvalue)
		tissue_vec <- c(tissue_vec, rep(tissue_type, length(temp_pvalue)))

	}
	df <- data.frame(pvalue=-log10(pvalue_vec+1e-30), tissue=factor(tissue_vec,levels=tissue_types))

	p <- ggplot(data=df, aes(x=tissue, y=pvalue, fill=tissue)) +
		geom_boxplot() +
		figure_theme() +
		theme(legend.position="none") + 
		labs(y="-log10(pvalue) for (trait, gene, tissue)", fill="", x="") +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
	return(p)
}


make_pip_count_barplot_stratefied_by_tissue_type <- function(ce_df) {
	pip_count <- c()
	tissue_vec <- c()

	num_samples <- dim(ce_df)[1]

	tissue_types <- as.character(sort(unique(ce_df$tissue)))


	for (tissue_iter in 1:length(tissue_types)) {
		tissue_type <- tissue_types[tissue_iter]

		indices <- ce_df$tissue == tissue_type
		temp_pip <- ce_df$pip[indices]
		temp_pip <- temp_pip[!is.na(temp_pip)]

		pip_count <- c(pip_count, sum(temp_pip > .5))
		tissue_vec <- c(tissue_vec, tissue_type)

	}
	df <- data.frame(pip_count=pip_count, tissue=factor(tissue_vec,levels=tissue_types))
	p<-ggplot(data=df, aes(x=tissue, y=pip_count)) +
  			geom_bar(stat="identity") +
  			figure_theme() + 
  			labs(y="Num genes PIP > .5", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
}

make_nominal_pvalue_thresh_count_barplot_stratefied_by_tissue_type <- function(ce_df) {
	pvalue_count_vec <- c()
	tissue_vec <- c()
	thresholds_vec <- c()

	num_samples <- dim(ce_df)[1]

	tissue_types <- as.character(sort(unique(ce_df$tissue)))

	pvalue_thresholds <- c(1e-4, 1e-6, 1e-8)
	for (pvalue_thresh_iter in 1:length(pvalue_thresholds)) {
		pvalue_thresh <- pvalue_thresholds[pvalue_thresh_iter]
		for (tissue_iter in 1:length(tissue_types)) {
			tissue_type <- tissue_types[tissue_iter]

			indices <- ce_df$tissue == tissue_type
			temp_pvalue <- ce_df$nominal_pvalue[indices]
			temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]

			pvalue_count_vec <- c(pvalue_count_vec, sum(temp_pvalue <= pvalue_thresh))
			tissue_vec <- c(tissue_vec, tissue_type)
			thresholds_vec <- c(thresholds_vec, pvalue_thresh)
		}
	}

	df <- data.frame(gene_count=pvalue_count_vec, tissue=factor(tissue_vec,levels=tissue_types), pvalue_threshold=factor(thresholds_vec, pvalue_thresholds))
	
	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun")

	p<-ggplot(data=df, aes(x=tissue, y=gene_count, fill=pvalue_threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="right") +
  			labs(y="Number of significant genes", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))


}


make_multivariate_pvalue_thresh_count_barplot_stratefied_by_tissue_type <- function(ce_df) {
	pvalue_count_vec <- c()
	tissue_vec <- c()
	thresholds_vec <- c()

	num_samples <- dim(ce_df)[1]

	tissue_types <- as.character(sort(unique(ce_df$tissue)))

	pvalue_thresholds <- c(1e-2, 1e-4, 1e-6, 1e-8, 1e-50)
	for (pvalue_thresh_iter in 1:length(pvalue_thresholds)) {
		pvalue_thresh <- pvalue_thresholds[pvalue_thresh_iter]
		for (tissue_iter in 1:length(tissue_types)) {
			tissue_type <- tissue_types[tissue_iter]

			indices <- ce_df$tissue == tissue_type
			temp_pvalue <- ce_df$multivariate_pvalue[indices]
			temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]

			pvalue_count_vec <- c(pvalue_count_vec, sum(temp_pvalue <= pvalue_thresh))
			tissue_vec <- c(tissue_vec, tissue_type)
			thresholds_vec <- c(thresholds_vec, pvalue_thresh)
		}
	}
	df <- data.frame(gene_count=pvalue_count_vec, tissue=factor(tissue_vec,levels=tissue_types), pvalue_threshold=factor(thresholds_vec, pvalue_thresholds))
	p<-ggplot(data=df, aes(x=tissue, y=gene_count, fill=pvalue_threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="top") +
  			labs(y="Num genes pvalue <= X", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))


}

make_multivariate_pvalue_thresh_rate_barplot_stratefied_by_tissue_type <- function(ce_df) {
	pvalue_count_vec <- c()
	tissue_vec <- c()
	thresholds_vec <- c()

	num_samples <- dim(ce_df)[1]

	tissue_types <- as.character(sort(unique(ce_df$tissue)))

	pvalue_thresholds <- c(1e-2, 1e-4, 1e-6, 1e-8)
	for (pvalue_thresh_iter in 1:length(pvalue_thresholds)) {
		pvalue_thresh <- pvalue_thresholds[pvalue_thresh_iter]
		for (tissue_iter in 1:length(tissue_types)) {
			tissue_type <- tissue_types[tissue_iter]

			indices <- ce_df$tissue == tissue_type
			temp_pvalue <- ce_df$multivariate_pvalue[indices]
			temp_pvalue <- temp_pvalue[!is.na(temp_pvalue)]

			pvalue_count_vec <- c(pvalue_count_vec, sum(temp_pvalue <= pvalue_thresh)/length(temp_pvalue))
			tissue_vec <- c(tissue_vec, tissue_type)
			thresholds_vec <- c(thresholds_vec, pvalue_thresh)
		}
	}
	df <- data.frame(gene_count=pvalue_count_vec, tissue=factor(tissue_vec,levels=tissue_types), pvalue_threshold=factor(thresholds_vec, pvalue_thresholds))
	p<-ggplot(data=df, aes(x=tissue, y=gene_count, fill=pvalue_threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="top") +
  			labs(y="Num genes pvalue <= X /\nNum tested genes", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))


}

make_pip_count_rate_barplot_stratefied_by_tissue_type <- function(ce_df) {
	pip_count <- c()
	tissue_vec <- c()

	num_samples <- dim(ce_df)[1]

	tissue_types <- as.character(sort(unique(ce_df$tissue)))


	for (tissue_iter in 1:length(tissue_types)) {
		tissue_type <- tissue_types[tissue_iter]

		indices <- ce_df$tissue == tissue_type
		temp_pip <- ce_df$pip[indices]
		temp_pip <- temp_pip[!is.na(temp_pip)]

		pip_count <- c(pip_count, sum(temp_pip > .5)/length(temp_pip))
		tissue_vec <- c(tissue_vec, tissue_type)

	}
	df <- data.frame(pip_count=pip_count, tissue=factor(tissue_vec,levels=tissue_types))
	p<-ggplot(data=df, aes(x=tissue, y=pip_count)) +
  			geom_bar(stat="identity") +
  			figure_theme() + 
  			labs(y="Num genes PIP > .5 /\n Num genes tested", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
}

make_fraction_of_time_each_tissue_is_causal_bar_blot_stratefied_by_tissue_type <- function(ce_df) {
	tissue_types <- as.character(sort(unique(ce_df$tissue)))
	genes <- as.character(sort(unique(ce_df$gene_name)))

	tissue_counts_hash <- hash()
	tissue_hits_hash <- hash()
	for (tissue_iter in 1:length(tissue_types)) {
		tissue_name <- tissue_types[tissue_iter]
		tissue_counts_hash[[tissue_name]] = 0
		tissue_hits_hash[[tissue_name]] = 0
	}


	for (gene_iter in 1:length(genes)) {
		current_gene <- genes[gene_iter]
		gene_indices <- ce_df$gene_name == current_gene

		gene_ce_df <- ce_df[gene_indices,]

		if (sum(!is.na(gene_ce_df$pip)) > 0) {

			if (sum(gene_ce_df$pip > .5) > 0) {
				observed_tissues <- as.character(gene_ce_df$tissue)
				observed_pips <- gene_ce_df$pip
				num_samp <- length(observed_tissues)


				for (samp_num in 1:num_samp) {
					observed_tissue <- observed_tissues[samp_num]
					observed_pip <- observed_pips[samp_num]
					tissue_counts_hash[[observed_tissue]] = tissue_counts_hash[[observed_tissue]] + 1
					if (observed_pip > .5) {
						tissue_hits_hash[[observed_tissue]] = tissue_hits_hash[[observed_tissue]] + 1
					}
				}
			}
		}


	}
	for (tissue_iter in 1:length(tissue_types)) {
		tissue_name <- tissue_types[tissue_iter]
		print(tissue_name)
		print(tissue_hits_hash[[tissue_name]])
		print(tissue_counts_hash[[tissue_name]])
	}
}


make_roc_curve <- function(ce_df, tissue_type) {
	valid_indices <- !is.na(ce_df$multivariate_pvalue)
	ce_df <- ce_df[valid_indices,]
	nominal_log_pvalue <- -log10(ce_df$nominal_pvalue + 1e-100)
	multivariate_log_pvalue <- -log10(ce_df$multivariate_pvalue + 1e-100)

	pos_labels = as.character(ce_df$tissue) == tissue_type
	neg_labels = as.character(ce_df$tissue) != tissue_type


	nominal_roc_obj2 = roc(response=pos_labels, predictor=nominal_log_pvalue, direction = ">")
	multivariate_roc_obj2 = roc(response=pos_labels, predictor=multivariate_log_pvalue, direction = ">")

	print(summary(nominal_roc_obj2))
	fpr <- c(1.0 - nominal_roc_obj2$specificities, 1.0 - multivariate_roc_obj2$specificities)
	tpr <- c(nominal_roc_obj2$sensitivities, multivariate_roc_obj2$sensitivities)
	model_names <- c(rep("nominal", length(nominal_roc_obj2$sensitivities)), rep("multivariate", length(multivariate_roc_obj2$sensitivities)))
	df <- data.frame(tpr=tpr, fpr=fpr, model_name=factor(model_names))

  	plotter <- ggplot(data=df, aes(x=fpr, y=tpr, group=model_name)) + geom_line(aes(colour=model_name)) + 
                labs(x="False positive rate", y="True positive rate",group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_abline() +
                figure_theme() 


	return(plotter)


}

make_roc_curve2 <- function(ce_df, tissue_type) {
	valid_indices <- !is.na(ce_df$multivariate_pvalue)
	ce_df <- ce_df[valid_indices,]



	#nom_sig_genes <- as.character(unique(ce_df$gene_name[as.character(ce_df$tissue) == tissue_type]))


	#valid_indices2 <- ce_df$gene_name %in% nom_sig_genes

	#ce_df <- ce_df[valid_indices2,]


	nom_sig_genes2 <- as.character(unique(ce_df$gene_name[ce_df$nominal_pvalue < 1e-5]))


	valid_indices3 <- ce_df$gene_name %in% nom_sig_genes2

	ce_df <- ce_df[valid_indices3,]	


	nominal_log_pvalue <- -log10(ce_df$nominal_pvalue + 1e-300)
	multivariate_log_pvalue <- -log10(ce_df$multivariate_pvalue + 1e-300)

	pos_labels = as.character(ce_df$tissue) == tissue_type
	neg_labels = as.character(ce_df$tissue) != tissue_type


	nominal_roc_obj2 = roc(response=pos_labels, predictor=-ce_df$nominal_pvalue, direction="<")
	multivariate_roc_obj2 = roc(response=pos_labels, predictor=-ce_df$multivariate_pvalue, direction="<")


	print(summary(nominal_log_pvalue[pos_labels]))
	print(summary(nominal_log_pvalue[neg_labels]))

	print(summary(multivariate_log_pvalue[pos_labels]))
	print(summary(multivariate_log_pvalue[neg_labels]))


	print(roc.test(multivariate_roc_obj2, nominal_roc_obj2)$p.value)

	fpr <- c(1.0 - nominal_roc_obj2$specificities, 1.0 - multivariate_roc_obj2$specificities)
	tpr <- c(nominal_roc_obj2$sensitivities, multivariate_roc_obj2$sensitivities)
	model_names <- c(rep("nominal", length(nominal_roc_obj2$sensitivities)), rep("multivariate", length(multivariate_roc_obj2$sensitivities)))
	df <- data.frame(tpr=tpr, fpr=fpr, model_name=factor(model_names))

  	plotter <- ggplot(data=df, aes(x=fpr, y=tpr, group=model_name)) + geom_line(aes(colour=model_name)) + 
                labs(x="False positive rate", y="True positive rate",group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_abline() +
                figure_theme() 


	return(plotter)


}

trait_name = args[1]
input_dir = args[2]
viz_dir = args[3]


# Load in data
causal_effects_file <- paste0(input_dir, "causal_effect_regression_merged.txt")
ce_df <- read.table(causal_effects_file, header=TRUE, sep="\t")

# Limit to data with a susie component 
susie_component_ce_df <- ce_df[(ce_df$num_trait_components != 0) & (ce_df$num_tissue_components != 0),]

# Load in data
twas_file <- paste0(input_dir, "UKB_460K.blood_WHITE_COUNT_marginal_twas_pvalues.txt")
twas_df <- read.table(twas_file, header=TRUE, sep="\t")



if (FALSE) {
# Make boxplot comparing distribution of -log10 pvalues for samples with susie component vs samples w/o a susie component
output_file <- paste0(viz_dir, "nominal_log10_pvalue_boxplot_stratefied_by_susie_component_status.pdf")
boxplot <- make_boxplot_comparing_nominal_log10_pvalues_depending_on_whether_there_was_a_susie_component(ce_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")


# Make boxplot comparing distribution of nominal -log10 pvalues for samples w/ susie component across tissues
output_file <- paste0(viz_dir, "nominal_log10_pvalue_boxplot_stratefied_by_tissue_type.pdf")
boxplot <- make_boxplot_comparing_nominal_log10_pvalues_stratefied_by_tissue_type(susie_component_ce_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")
}


# Make baplot showing number samples with nominal pvalue < thresh in each tissue
output_file <- paste0(viz_dir, "twas_nominal_pvalue_thresh_count_bar_blot_stratefied_by_tissue_type.pdf")
boxplot <- make_nominal_pvalue_thresh_count_barplot_stratefied_by_tissue_type(twas_df)
ggsave(boxplot, file=output_file, width=7.9, height=4.3, units="in")

print("DONE")

# Make ROC CURVE
tissue_type <- "Whole_Blood"
output_file <- paste0(viz_dir, "roc_curve_", tissue_type, ".pdf")
#roc <- make_roc_curve(susie_component_ce_df, tissue_type)
#ggsave(roc, file=output_file, width=7.2, height=5.0, units="in")

# Make ROC CURVE
tissue_type <- "Whole_Blood"
output_file <- paste0(viz_dir, "roc_curve2_", tissue_type, ".pdf")
roc <- make_roc_curve2(susie_component_ce_df, tissue_type)
ggsave(roc, file=output_file, width=7.2, height=5.0, units="in")


# Make baplot showing fraction of times each tissue was causal for an observed causal variant 
output_file <- paste0(viz_dir, "fraction_of_time_each_tissue_is_causal_bar_blot_stratefied_by_tissue_type.pdf")
boxplot <- make_fraction_of_time_each_tissue_is_causal_bar_blot_stratefied_by_tissue_type(susie_component_ce_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")


# Make baplot showing number samples with multivariate pvalue < thresh in each tissue
output_file <- paste0(viz_dir, "multivariate_pvalue_thresh_count_bar_blot_stratefied_by_tissue_type.pdf")
boxplot <- make_multivariate_pvalue_thresh_count_barplot_stratefied_by_tissue_type(susie_component_ce_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")


# Make baplot showing number samples with rate of multivariate pvalue < thresh in each tissue
output_file <- paste0(viz_dir, "multivariate_pvalue_thresh_rate_bar_blot_stratefied_by_tissue_type.pdf")
boxplot <- make_multivariate_pvalue_thresh_rate_barplot_stratefied_by_tissue_type(susie_component_ce_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")



# Make baplot showing number samples with nominal pvalue < thresh in each tissue
output_file <- paste0(viz_dir, "nominal_pvalue_thresh_count_bar_blot_stratefied_by_tissue_type.pdf")
boxplot <- make_nominal_pvalue_thresh_count_barplot_stratefied_by_tissue_type(susie_component_ce_df)
ggsave(boxplot, file=output_file, width=7.9, height=5.0, units="in")




# Make baplot showing number samples with PIP > Thresh in each tissue
output_file <- paste0(viz_dir, "pip_count_bar_blot_stratefied_by_tissue_type.pdf")
boxplot <- make_pip_count_barplot_stratefied_by_tissue_type(susie_component_ce_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")

# Make baplot showing rate of number samples with PIP > Thresh in each tissue
output_file <- paste0(viz_dir, "pip_count_rate_bar_blot_stratefied_by_tissue_type.pdf")
boxplot <- make_pip_count_rate_barplot_stratefied_by_tissue_type(susie_component_ce_df)
ggsave(boxplot, file=output_file, width=7.2, height=5.0, units="in")


