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

extract_fraction_of_twas_nominal_associations <- function(gtex_fusion_associations_dir, tissue_names, trait_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	thresh_arr <- c()

	thresholds <- c(1e-4, 1e-6, 1e-8, 1e-10)


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			pvalz <- c()
			for (chrom_num in 1:22) {

				twas_assoc_file <- paste0(gtex_fusion_associations_dir, tissue_name, ".", trait_name, "_", chrom_num, ".dat")
				twas_df <- read.table(twas_assoc_file, header=TRUE)
				pvalz <- c(pvalz, twas_df$TWAS.P)
			}

			for (threshold_iter in 1:length(thresholds)) {
				threshold <- thresholds[threshold_iter]
				counts = sum(pvalz <= threshold)/length(pvalz)
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				thresh_arr <- c(thresh_arr, threshold)
				count_arr <- c(count_arr, counts)
			}
		}

	}

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), pvalue_threshold=factor(thresh_arr), fraction_nominal_twas_associations=count_arr)
}

extract_fraction_of_twas_multivariate_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	thresh_arr <- c()

	thresholds <- c(.5, .7, .9)


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		file_name <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

		df <- read.table(file_name, header=TRUE)
		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- as.character(tissue_names[tissue_iter]) 

			tissue_df <- df[df$tissue==tissue_name,]
			for (threshold_iter in 1:length(thresholds)) {
				threshold <- thresholds[threshold_iter]

				counts <- sum(tissue_df$multivariate_twas_pip >= threshold)/length(tissue_df$multivariate_twas_pip)
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				thresh_arr <- c(thresh_arr, threshold)
				count_arr <- c(count_arr, counts)
			}
		}
	}
	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), pvalue_threshold=factor(thresh_arr), fraction_nominal_twas_associations=count_arr)
	return(df)
}

extract_number_of_twas_multivariate_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	thresh_arr <- c()

	thresholds <- c(.5, .7, .9)


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		file_name <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

		df <- read.table(file_name, header=TRUE)
		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- as.character(tissue_names[tissue_iter]) 

			tissue_df <- df[df$tissue==tissue_name,]
			for (threshold_iter in 1:length(thresholds)) {
				threshold <- thresholds[threshold_iter]

				counts <- sum(tissue_df$multivariate_twas_pip >= threshold)
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				thresh_arr <- c(thresh_arr, threshold)
				count_arr <- c(count_arr, counts)
			}
		}
	}
	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), pvalue_threshold=factor(thresh_arr), number_nominal_twas_associations=count_arr)
	return(df)
}
extract_number_of_max_twas_multivariate_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	thresh_arr <- c()

	thresholds <- c(.5, .7, .9)


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		file_name <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

		df <- read.table(file_name, header=TRUE)
		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- as.character(tissue_names[tissue_iter]) 

			tissue_df <- df[df$tissue==tissue_name,]

			trait_components = unique(tissue_df$trait_component)
			max_pips <- c()

			for (trait_component_iter in 1:length(trait_components)) {
				trait_component <- trait_components[trait_component_iter]

				max_pip = max(tissue_df$multivariate_twas_pip[tissue_df$trait_component==trait_component])
				max_pips <- c(max_pips, max_pip)
			}


			for (threshold_iter in 1:length(thresholds)) {
				threshold <- thresholds[threshold_iter]

				counts <- sum(max_pips >= threshold)/length(max_pips)
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				thresh_arr <- c(thresh_arr, threshold)
				count_arr <- c(count_arr, counts)
			}
		}
	}
	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), pvalue_threshold=factor(thresh_arr), number_nominal_twas_associations=count_arr)
	return(df)
}


extract_number_of_competitive_max_twas_nominal_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)
	tissue_names = as.character(tissue_names)

	pseudocount = .005


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(trait_name)

		file_name <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

		df <- read.table(file_name, header=TRUE)

		tissue_counts <- rep(0, num_tiss)
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- tissue_names[tissue_iter]
			num_genes <- sum(df$tissue==tissue_name)
			tissue_counts[tissue_iter] = num_genes
		}
		tissue_counts = tissue_counts/mean(tissue_counts)
		#print(tissue_names)
		#print(tissue_counts)


		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name
				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(min(component_df$nominal_twas_pvalue[indices]) < 1e-6)
				}
			
			}

			tissues_un_normalized = (tissues_un_normalized) + pseudocount
			tissues_normalized = tissues_un_normalized/sum(tissues_un_normalized)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	print((df))
	return(df)
}


extract_number_of_genes_per_trait_component <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	count_arr <- c()
	comp_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(trait_name)

		file_name <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

		df <- read.table(file_name, header=TRUE)
		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]

			trait_arr <- c(trait_arr, trait_name)
			comp_arr <- c(comp_arr, trait_component)
			count_arr <- c(count_arr, dim(component_df)[1])

		}
	}
	df <- data.frame(trait=factor(trait_arr), component=factor(comp_arr), num_genes=count_arr)
		p<-ggplot(df, aes(x=trait, y=num_genes, fill=trait)) +
  			geom_violin(trim=FALSE) + figure_theme() +
  			theme(legend.position="none")

  		return(p)

 }


extract_number_of_competitive_max_twas_multivariate_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)
	tissue_names = as.character(tissue_names)

	pseudocount = .005


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		print(trait_name)

		file_name <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

		df <- read.table(file_name, header=TRUE)

		tissue_counts <- rep(0, num_tiss)
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- tissue_names[tissue_iter]
			num_genes <- sum(df$tissue==tissue_name)
			tissue_counts[tissue_iter] = num_genes
		}
		tissue_counts = tissue_counts/mean(tissue_counts)


		
		trait_components = as.character(unique(df$trait_component))
		for (trait_component_iter in 1:length(trait_components)) {
			trait_component <- trait_components[trait_component_iter]
			component_df <- df[df$trait_component==trait_component,]
			# initialize
			tissues_un_normalized <- rep(0.0, num_tiss)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]

				indices <- component_df$tissue==tissue_name

				if (sum(indices) > 0) {
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(component_df$multivariate_twas_pip[indices]) > .5)
				}
			}
			tissues_un_normalized = (tissues_un_normalized) + pseudocount

			tissues_normalized = tissues_un_normalized/sum(tissues_un_normalized)
			for (tissue_iter in 1:num_tiss) {
				tissue_name <- tissue_names[tissue_iter]
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				count_arr <- c(count_arr, tissues_normalized[tissue_iter])
				comp_arr <- c(comp_arr, trait_component)
			}			

		}
	}

	df <- data.frame(trait=trait_arr, tissue=(tissue_arr), trait_component=comp_arr, number_nominal_twas_associations=count_arr)

	trait_arr <- c()
	tissue_arr <- c()
	fraction_arr <- c()

	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

		df_trait <- df[df$trait==trait_name,]
		for (tissue_iter in 1:num_tiss) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			probs <- df_trait$number_nominal_twas_associations[df_trait$tissue==tissue_name]
			frac = sum(probs)/length(probs)

			trait_arr <- c(trait_arr, trait_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			fraction_arr <- c(fraction_arr, frac)

		}
	}	

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), fraction_nominal_twas_associations=fraction_arr)

	print((df))
	return(df)
}






extract_number_of_twas_nominal_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	thresh_arr <- c()

	thresholds <- c(1e-4, 1e-6, 1e-8, 1e-10)


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		file_name <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_component_organized_multivariate_twas_overlaps.txt")

		df <- read.table(file_name, header=TRUE)
		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- as.character(tissue_names[tissue_iter]) 

			tissue_df <- df[df$tissue==tissue_name,]
			for (threshold_iter in 1:length(thresholds)) {
				threshold <- thresholds[threshold_iter]

				counts <- sum(tissue_df$nominal_twas_pvalue <=threshold)
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				thresh_arr <- c(thresh_arr, threshold)
				count_arr <- c(count_arr, counts)
			}
		}
	}
	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), pvalue_threshold=factor(thresh_arr), number_nominal_twas_associations=count_arr)
	return(df)
}



extract_number_of_twas_nominal_associations <- function(gtex_fusion_associations_dir, tissue_names, trait_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	thresh_arr <- c()

	thresholds <- c(1e-4, 1e-6, 1e-8, 1e-10)


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])
		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- as.character(tissue_names[tissue_iter])
			pvalz <- c()
			for (chrom_num in 1:22) {

				twas_assoc_file <- paste0(gtex_fusion_associations_dir, tissue_name, ".", trait_name, "_", chrom_num, ".dat")
				twas_df <- read.table(twas_assoc_file, header=TRUE)
				pvalz <- c(pvalz, twas_df$TWAS.P)
			}

			for (threshold_iter in 1:length(thresholds)) {
				threshold <- thresholds[threshold_iter]
				counts = sum(pvalz <= threshold)
				trait_arr <- c(trait_arr, trait_name)
				tissue_arr <- c(tissue_arr, tissue_name)
				thresh_arr <- c(thresh_arr, threshold)
				count_arr <- c(count_arr, counts)
			}
		}

	}

	df <- data.frame(trait=trait_arr, tissue=factor(tissue_arr), pvalue_threshold=factor(thresh_arr), number_nominal_twas_associations=count_arr)
}

make_nominal_pip_thresh_bar_plot <- function(df, trait_name) {
	df = df[df$trait==trait_name, ]

	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn")

	df_small = df[df$pvalue_threshold=="0.5",]

	indices <- order(df_small$number_nominal_twas_associations)

	tissue_order = as.character(df_small$tissue)[indices]

	df$tissue = factor(df$tissue, levels=tissue_order)

	p<-ggplot(data=df, aes(x=tissue, y=number_nominal_twas_associations, fill=pvalue_threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="right") +
  			labs(y="Number of significant genes", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10), axis.text.y=element_text(size=10))

  	return(p)

}

make_nominal_pvalue_thresh_bar_plot <- function(df, trait_name) {
	df = df[df$trait==trait_name, ]

	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn")

	df_small = df[df$pvalue_threshold=="1e-10",]

	indices <- order(df_small$number_nominal_twas_associations)

	tissue_order = as.character(df_small$tissue)[indices]

	df$tissue = factor(df$tissue, levels=tissue_order)

	p<-ggplot(data=df, aes(x=tissue, y=number_nominal_twas_associations, fill=pvalue_threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="right") +
  			labs(y="Number of significant genes", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10), axis.text.y=element_text(size=10))

  	return(p)

}

make_nominal_pvalue_thresh_fraction_bar_plot <- function(df, trait_name) {
	df = df[df$trait==trait_name, ]

	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn")

	df_small = df[df$pvalue_threshold=="1e-10",]

	indices <- order(df_small$fraction_nominal_twas_associations)

	tissue_order = as.character(df_small$tissue)[indices]

	df$tissue = factor(df$tissue, levels=tissue_order)



	p<-ggplot(data=df, aes(x=tissue, y=fraction_nominal_twas_associations, fill=pvalue_threshold)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			theme(legend.position="right") +
  			labs(y="Fraction of heritible genes significant", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))

  	return(p)

}


make_comp_nominal_pvalue_thresh_fraction_bar_plot <- function(df, trait_name) {
	df = df[df$trait==trait_name, ]

	df$tissue <- recode(df$tissue, Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn")

	indices <- order(df$fraction_nominal_twas_associations)

	tissue_order = as.character(df$tissue)[indices]

	df$tissue = factor(df$tissue, levels=tissue_order)



	p<-ggplot(data=df, aes(x=tissue, y=fraction_nominal_twas_associations)) +
  			geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
  			figure_theme() + 
  			labs(y="Fraction of mediated trait components", x="") +
  			theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))

  	return(p)

}



gtex_tissue_file <- args[1]
gtex_fusion_associations_dir <- args[2]
gtex_fusion_multivariate_associations_dir <- args[3]
visualize_twas_results_dir <- args[4]



# Load in gtex tissues
tissue_df <- read.table(gtex_tissue_file,header=TRUE)
tissue_names <- tissue_df$tissue_name

# names of traits 

trait_names <- c("blood_WHITE_COUNT", "lung_FEV1FVCzSMOKE", "body_WHRadjBMIz", "bp_DIASTOLICadjMEDz")


competitive_max_twas_nominal_associations_fraction_nearby_trait_components_df <- extract_number_of_competitive_max_twas_nominal_associations_nearby_trait_components(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names)

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_competitive_twas_nominal_pvalue_thresh_fraction_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(competitive_max_twas_nominal_associations_fraction_nearby_trait_components_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}


num_competitive_max_twas_multivariate_associations_nearby_trait_components_df <- extract_number_of_competitive_max_twas_multivariate_associations_nearby_trait_components(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names)
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_competitive_twas_multivariate_pvalue_thresh_fraction_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(num_competitive_max_twas_multivariate_associations_nearby_trait_components_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}

print("DONE")


if (FALSE) {
num_genes_per_trait_component_violinplot <- extract_number_of_genes_per_trait_component(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names)
ggsave(num_genes_per_trait_component_violinplot, file=paste0(visualize_twas_results_dir, "num_genes_per_component_violin.pdf"), width=7.9, height=4.3, units="in")
}



print("DONE")



# Extract number of twas nominal associations in each trait, tissue at various p-value thresholds nearby trait components
#num_twas_nominal_associations_nearby_trait_components_df <- extract_number_of_twas_nominal_associations_nearby_trait_components(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names)
num_twas_multivariate_associations_nearby_trait_components_df <- extract_number_of_twas_multivariate_associations_nearby_trait_components(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names)
#fraction_twas_multivariate_associations_nearby_trait_components_df <- extract_fraction_of_twas_multivariate_associations_nearby_trait_components(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names)
#num_max_twas_multivariate_associations_nearby_trait_components_df <- extract_number_of_max_twas_multivariate_associations_nearby_trait_components(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names)



for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_twas_multivariate_pip_thresh_count_nearby_trait_components_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_nominal_pip_thresh_bar_plot(num_twas_multivariate_associations_nearby_trait_components_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}

print("DONE")

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_twas_nominal_pvalue_thresh_count_nearby_trait_components_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_nominal_pvalue_thresh_bar_plot(num_twas_nominal_associations_nearby_trait_components_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_twas_multivariate_max_pip_thresh_count_nearby_trait_components_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_nominal_pvalue_thresh_bar_plot(num_max_twas_multivariate_associations_nearby_trait_components_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_twas_multivariate_pip_thresh_fraction_nearby_trait_components_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_nominal_pvalue_thresh_fraction_bar_plot(fraction_twas_multivariate_associations_nearby_trait_components_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}



# Extract number of twas nominal associations in each trait, tissue at various p-value thresholds
num_twas_nominal_associations_df <- extract_number_of_twas_nominal_associations(gtex_fusion_associations_dir, tissue_names, trait_names)
fraction_twas_nominal_associations_df <- extract_fraction_of_twas_nominal_associations(gtex_fusion_associations_dir, tissue_names, trait_names)


for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_twas_nominal_pvalue_thresh_count_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_nominal_pvalue_thresh_bar_plot(num_twas_nominal_associations_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}

for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	output_file <- paste0(visualize_twas_results_dir, trait_name, "_twas_nominal_pvalue_thresh_fraction_bar_blot_stratefied_by_tissue_type.pdf")
	barplot <- make_nominal_pvalue_thresh_fraction_bar_plot(fraction_twas_nominal_associations_df, trait_name)
	ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")
}



