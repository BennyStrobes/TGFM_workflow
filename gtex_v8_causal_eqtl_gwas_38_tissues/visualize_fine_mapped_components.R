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



get_causal_tissues_corresponding_to_trait <- function(trait_name) {
	if (trait_name == "blood_WHITE_COUNT") {
		causal_tissue =c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "body_WHRadjBMIz") {
		causal_tissue = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum")
		single_causal_tissue = "Adipose_Subcutaneous"
	}
	if (trait_name == "biochemistry_Cholesterol") {
		causal_tissue = c("Liver")
		single_causal_tissue = "Liver"
	}
	if (trait_name == "bp_DIASTOLICadjMEDz") {
		causal_tissue = c("Artery_Heart", "Artery_Tibial", "Heart_Atrial_Appendage", "Heart_Left_Ventricle")
		single_causal_tissue = "Artery_Tibial"
	}
	if (trait_name == "blood_EOSINOPHIL_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "blood_RBC_DISTRIB_WIDTH") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "blood_RED_COUNT") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Whole_Blood"
	}
	if (trait_name == "cov_EDU_COLLEGE") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
		single_causal_tissue = "Brain_Cerebellum"
	}
	if (trait_name == "mental_NEUROTICISM") {
		causal_tissue <- c("Brain_BasalGanglia", "Brain_Cerebellum", "Brain_Cortex", "Brain_Limbic", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra")
		single_causal_tissue = "Brain_Cerebellum"
	}
	if (trait_name == "disease_ALLERGY_ECZEMA_DIAGNOSED") {
		causal_tissue <- c("Whole_Blood", "Cells_EBV_transformed_lymphocytes")
		single_causal_tissue = "Cells_EBV_transformed_lymphocytes"
	}

	ret <- list()
	ret[[1]] = causal_tissue
	ret[[2]] = single_causal_tissue

	return(ret)
}


fine_mapped_manhatten_plot_orig <- function(df) {
	
	df$log_10_pvalue = -log10(df$p_value)
	df$position = df$position/1000000.0

	p <- ggplot(df, aes(x=position, y=log_10_pvalue,  color=fm_annotation)) +
  		geom_point() +
  		figure_theme() + 
  		theme(legend.position="bottom") + 
  		labs(color="", y="-log10(p-value)", x="Position (MB)") +
  		guides(color = guide_legend( nrow=3, byrow=TRUE)) +
    	theme(legend.key.size = unit(.8, 'cm'), legend.title = element_text(size=9),legend.text = element_text(size=9))

  	return(p)
}



fine_mapped_manhatten_plot <- function(df) {
	
	df$log_10_pvalue = -log10(df$p_value)
	df$position = df$position/1000000.0


	df_null = df[is.na(df$fm_annotation),]
	df_comp = df[!is.na(df$fm_annotation),]


	p <- ggplot() +
  		geom_point(data=df_null, aes(x=position, y=log_10_pvalue), color='grey', size=.5) +
    	geom_point(data=df_comp, aes(x=position, y=log_10_pvalue, color=fm_annotation), size=2.6) +
  		figure_theme() + 
  		theme(legend.position="top") + 
  		labs(color="", y="-log10(p-value)", x="Position (MB)") +
  		guides(color = guide_legend( nrow=2, byrow=TRUE)) +
    	theme(legend.key.size = unit(.8, 'cm'), legend.title = element_text(size=9),legend.text = element_text(size=9))

  	return(p)
}



component_contains_causal_annotation <- function(component_annotations, single_causal_tissue) {
	pass_bool = FALSE
	for (anno_iter in 1:length(component_annotations)) {
		annotation <- component_annotations[anno_iter]
		if (grepl(single_causal_tissue,annotation)) {
			pass_bool = TRUE
		}

	}
	return(pass_bool)
}

component_contains_two_annotations <- function(component_annotations, single_causal_tissue1, single_causal_tissue2) {
	pass_bool1 = FALSE
	pass_bool2 = FALSE
	for (anno_iter in 1:length(component_annotations)) {
		annotation <- component_annotations[anno_iter]
		if (grepl(single_causal_tissue1,annotation)) {
			pass_bool1 = TRUE
		}
		if (grepl(single_causal_tissue2,annotation)) {
			pass_bool2 = TRUE
		}

	}
	pass_bool=FALSE
	if (pass_bool1) {
		if (pass_bool2) {
			pass_bool = TRUE
		}
	}

	return(pass_bool)
}





trait_file <- args[1]
gtex_pseudotissue_file <- args[2]
pseudotissue_gtex_rss_multivariate_twas_dir <- args[3]
output_dir <- args[4]
gtex_tissue_colors_file <- args[5]


# Model parameters
fusion_weights="False"
gene_version = "cis_heritable_genes"

# Name of trait
trait_name <- "blood_WHITE_COUNT"
trait_name <- "bp_DIASTOLICadjMEDz"



# Get causal tissue and single_causal tissue corresponding to this trait
causal_tissue_list <- get_causal_tissues_corresponding_to_trait(trait_name)
causal_tissue <- causal_tissue_list[[1]]
single_causal_tissue <- causal_tissue_list[[2]]

# Load in organized TGFM results
trait_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_component_organized_joint_tgfm_results.txt")
trait_df <- read.table(trait_file, header=TRUE)


# Get unique components for this trait
components <- unique(as.character(trait_df$trait_component))

for (component_iter in 1:length(components)) {
	component_name <- components[component_iter]
	namer <- paste0(trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_", component_name)
	component_manhatten_info_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_dir, trait_name, "_", gene_version, "_fusion_weights_", fusion_weights, "_", component_name, "_manhatten_plot_info.txt")
	df <- read.table(component_manhatten_info_file, header=TRUE)

	component_annotations <- unique(as.character(df$fm_annotation))
	if (component_contains_causal_annotation(component_annotations, single_causal_tissue)) {
		fm_manhatten_plot <- fine_mapped_manhatten_plot(df)
		output_file <- paste0(output_dir, namer, "_fm_manhatten_plot.pdf")
		ggsave(fm_manhatten_plot, file=output_file, width=7.2, height=4.0, units="in")

	}
}










