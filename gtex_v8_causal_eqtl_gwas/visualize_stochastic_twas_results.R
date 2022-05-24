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


nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot <- function(stochastic_twas_df, fusion_twas_df) {
	df = data.frame(fusion_twas_z=fusion_twas_df$nominal_twas_zscore, susie_pmces_twas_z=stochastic_twas_df$nominal_twas_zscore)
	corry <- cor(df$fusion_twas_z, df$susie_pmces_twas_z)
	p <- ggplot(df, aes(x=fusion_twas_z, y=susie_pmces_twas_z)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline() +
		labs(title=paste0("Adipose Subcutaneous tissue / pearson correlation: ", corry))
	return(p)	
}

abs_diff_nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot_by_stochastic_z_score_sd <- function(stochastic_twas_df, fusion_twas_df) {
	df = data.frame(absolute_twas_z_diff=abs(fusion_twas_df$nominal_twas_zscore - stochastic_twas_df$nominal_twas_zscore), stochastic_twas_z_sd=stochastic_twas_df$nominal_stochastic_zscore_sd)
	p <- ggplot(df, aes(x=stochastic_twas_z_sd, y=absolute_twas_z_diff)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline() +
		labs(title=paste0("Adipose Subcutaneous tissue"))
	return(p)	

}


nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot_colored_by_stochastic_z_score_sd <- function(stochastic_twas_df, fusion_twas_df) {

	df = data.frame(fusion_twas_z=fusion_twas_df$nominal_twas_zscore, susie_pmces_twas_z=stochastic_twas_df$nominal_twas_zscore, stochastic_twas_z_sd=log(1+stochastic_twas_df$nominal_stochastic_zscore_sd))
	p <- ggplot(df, aes(x=fusion_twas_z, y=susie_pmces_twas_z, color=stochastic_twas_z_sd)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline()
	return(p)	
}


make_nominal_zscore_stochastic_mean_zscore_scatterplot <- function(stochastic_twas_df) {
	p <- ggplot(stochastic_twas_df, aes(x=nominal_stochastic_zscore_mean, y=nominal_twas_zscore, color=nominal_stochastic_zscore_sd)) + geom_point(size=.1) +
		figure_theme() +
		geom_abline() + 
		labs(x="empirical_mean_susie_twas_z", y="susie_pmces_twas_z", color="empirical_sd_susie_twas_z")
	return(p)
}

make_nominal_zscore_stochastic_mean_zscore_scatterplot_colored_by_vmr <- function(stochastic_twas_df) {
	stochastic_twas_df$vmr = ((stochastic_twas_df$nominal_stochastic_zscore_sd)^2)/(stochastic_twas_df$multivariate_twas_stochastically_averaged_pip)
	
	stochastic_twas_df$vmr[stochastic_twas_df$vmr > 100] = 100
	stochastic_twas_df$vmr = stochastic_twas_df$vmr
	p <- ggplot(stochastic_twas_df, aes(x=abs(nominal_stochastic_zscore_mean-nominal_twas_zscore), y=vmr)) + geom_point(size=.001) +
		figure_theme() +
		geom_abline() + 
		labs(x="abs(empirical_mean_susie_twas_z-empirical_vmr_susie_twas_z)", y="empirical_vmr_susie_twas_z")
	return(p)
}

make_nominal_pip_stochastic_mean_pip_scatterplot <-function(stochastic_twas_df) {
	print(cor(stochastic_twas_df$multivariate_twas_stochastically_averaged_pip, stochastic_twas_df$multivariate_twas_pip))
	p <- ggplot(stochastic_twas_df, aes(x=multivariate_twas_stochastically_averaged_pip, y=multivariate_twas_pip)) + geom_point() +
		figure_theme() +
		geom_abline()
	return(p)
}

make_nominal_pip_stochastic_mean_pip_scatterplot_colored_by_whole_blood <-function(stochastic_twas_df) {
	stochastic_twas_df$whole_blood_boolean = stochastic_twas_df$tissue=="Whole_Blood"
	p <- ggplot(stochastic_twas_df, aes(x=multivariate_twas_stochastically_averaged_pip, y=multivariate_twas_pip, color=whole_blood_boolean)) + geom_point() +
		figure_theme() +
		geom_abline()
	return(p)
}

make_nominal_pip_stochastic_mean_pip_scatterplot_colored_by_zscore_sd <-function(stochastic_twas_df) {
	p <- ggplot(stochastic_twas_df, aes(x=multivariate_twas_stochastically_averaged_pip, y=multivariate_twas_pip, color=nominal_stochastic_zscore_sd)) + geom_point() +
		figure_theme() +
		geom_abline()
	return(p)
}

make_nominal_pip_stochastic_mean_pip_scatterplot_colored_by_vmr <-function(stochastic_twas_df) {
	stochastic_twas_df$vmr = ((stochastic_twas_df$nominal_stochastic_zscore_sd)^2)/(stochastic_twas_df$multivariate_twas_stochastically_averaged_pip)
	
	stochastic_twas_df$vmr[stochastic_twas_df$vmr > 100] = 100
	stochastic_twas_df$vmr = stochastic_twas_df$vmr

	corry = cor(stochastic_twas_df$multivariate_twas_stochastically_averaged_pip, stochastic_twas_df$multivariate_twas_pip)

	p <- ggplot(stochastic_twas_df, aes(x=multivariate_twas_stochastically_averaged_pip, y=multivariate_twas_pip, color=vmr)) + geom_point(size=.01) +
		figure_theme() +
		geom_abline() + 
		labs(color="Emperical pip variance/\nempirical pip mean", title=paste0("spearman correlation: ", corry))
	return(p)
}



extract_number_of_competitive_max_stochastic_twas_multivariate_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)
	tissue_names = as.character(tissue_names)

	pseudocount = .0001


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(max(component_df$multivariate_twas_stochastically_averaged_pip[indices]) > .5)
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + (sum(component_df$multivariate_twas_stochastically_averaged_pip[indices]))
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
	return(df)
}
extract_number_of_competitive_max_twas_nominal_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)
	tissue_names = as.character(tissue_names)

	pseudocount = .0001


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

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
					tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + 1.0*(min(component_df$nominal_twas_pvalue[indices]) < 1e-6)
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + (sum(component_df$multivariate_twas_pip[indices])/length(indices))
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
	return(df)
}


extract_number_of_competitive_max_twas_multivariate_associations_nearby_trait_components <- function(gtex_fusion_multivariate_associations_dir, trait_names, tissue_names) {
	trait_arr <- c()
	tissue_arr <- c()
	count_arr <- c()
	comp_arr <- c()


	num_tiss <- length(tissue_names)
	tissue_names = as.character(tissue_names)

	pseudocount = .0001


	for (trait_iter in 1:length(trait_names)) {
		trait_name <- as.character(trait_names[trait_iter])

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
					#tissues_un_normalized[tissue_iter] = tissues_un_normalized[tissue_iter] + (sum(component_df$multivariate_twas_pip[indices])/length(indices))
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

	return(df)
}

make_comp_nominal_pvalue_thresh_fraction_bar_plot <- function(df) {

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

make_comp_twas_model_comparison_boxplot <- function(competitive_df) {
	competitive_df$whole_blood_boolean = as.character(competitive_df$tissue)=="Whole_Blood"
	p <- ggplot(data=competitive_df, aes(x=model, y=fraction_nominal_twas_associations, fill=whole_blood_boolean)) +
		geom_boxplot() +
		figure_theme() +
		labs(y="Fraction of mediated trait components")

	return(p)
}



gtex_tissue_file <- args[1]
gtex_fusion_multivariate_stochastic_associations_dir <- args[2]
gtex_fusion_multivariate_associations_dir <- args[3]
visualize_twas_results_dir <- args[4]


# Load in gtex tissues
tissue_df <- read.table(gtex_tissue_file,header=TRUE)
tissue_names <- tissue_df$pseudotissue_name# names of traits 

# Names of traits
trait_names <- c("blood_WHITE_COUNT")

nominal_twas_competitive_df <- extract_number_of_competitive_max_twas_nominal_associations_nearby_trait_components(gtex_fusion_multivariate_stochastic_associations_dir, trait_names, tissue_names)
stochastic_multivariate_twas_competitive_df <- extract_number_of_competitive_max_stochastic_twas_multivariate_associations_nearby_trait_components(gtex_fusion_multivariate_stochastic_associations_dir, trait_names, tissue_names)
multivariate_twas_competitive_df <- extract_number_of_competitive_max_twas_multivariate_associations_nearby_trait_components(gtex_fusion_multivariate_stochastic_associations_dir, trait_names, tissue_names)

#competitive_df <- data.frame(model=c(rep("univariate_twas", 23), rep("multivariate_twas", 23), rep("multivariate_stochastic_twas"), 23),trait=c(nominal_twas_competitive_df$trait, stochastic_multivariate_twas_competitive_df$trait, multivariate_twas_competitive_df$trait),fraction_nominal_twas_associations=c(nominal_twas_competitive_df$fraction_nominal_twas_associations, stochastic_multivariate_twas_competitive_df$fraction_nominal_twas_associations, multivariate_twas_competitive_df$fraction_nominal_twas_associations), tissue=c(nominal_twas_competitive_df$tissue, stochastic_multivariate_twas_competitive_df$tissue, multivariate_twas_competitive_df$tissue))
competitive_df <- data.frame(model=as.character(c(rep("univariate_twas", 23), rep("multivariate_stochastic_twas", 23), rep("multivariate_twas", 23))), trait=as.character(c(as.character(nominal_twas_competitive_df$trait), as.character(stochastic_multivariate_twas_competitive_df$trait), as.character(multivariate_twas_competitive_df$trait))),fraction_nominal_twas_associations=c(nominal_twas_competitive_df$fraction_nominal_twas_associations, stochastic_multivariate_twas_competitive_df$fraction_nominal_twas_associations, multivariate_twas_competitive_df$fraction_nominal_twas_associations), tissue=as.character(c(as.character(nominal_twas_competitive_df$tissue), as.character(stochastic_multivariate_twas_competitive_df$tissue), as.character(multivariate_twas_competitive_df$tissue))))
competitive_df$model = factor(competitive_df$model, levels=c("univariate_twas", "multivariate_twas", "multivariate_stochastic_twas"))


output_file <- paste0(visualize_twas_results_dir, "competitive_twas_various_models.pdf")
barplot <- make_comp_twas_model_comparison_boxplot(competitive_df)
ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")




output_file <- paste0(visualize_twas_results_dir, "competitive_twas_multivariate_pvalue_thresh_fraction_bar_blot_stratefied_by_tissue_type.pdf")
barplot <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(multivariate_twas_competitive_df)
ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")


output_file <- paste0(visualize_twas_results_dir, "competitive_stochastic_twas_multivariate_pip_thresh_fraction_bar_blot_stratefied_by_tissue_type.pdf")
barplot <- make_comp_nominal_pvalue_thresh_fraction_bar_plot(stochastic_multivariate_twas_competitive_df)
ggsave(barplot, file=output_file, width=7.9, height=4.3, units="in")




# Extract twas stochastic df
stochastic_twas_df <- read.table(paste0(gtex_fusion_multivariate_stochastic_associations_dir, "blood_WHITE_COUNT_component_organized_multivariate_twas_overlaps.txt"), header=TRUE)
stochastic_twas_df$nominal_twas_zscore[stochastic_twas_df$nominal_twas_zscore >45]=45
stochastic_twas_df$nominal_twas_zscore[stochastic_twas_df$nominal_twas_zscore < -45]=-45
stochastic_twas_df$nominal_stochastic_zscore_mean[stochastic_twas_df$nominal_stochastic_zscore_mean >45]=45
stochastic_twas_df$nominal_stochastic_zscore_mean[stochastic_twas_df$nominal_stochastic_zscore_mean < -45]=-45
stochastic_twas_df$nominal_stochastic_zscore_sd[stochastic_twas_df$nominal_stochastic_zscore_sd > 35] = 35



fusion_twas_df <- read.table(paste0(gtex_fusion_multivariate_associations_dir, "blood_WHITE_COUNT_component_organized_multivariate_twas_overlaps.txt"), header=TRUE)





fusion_twas_df_adipose_sub <- fusion_twas_df[as.character(fusion_twas_df$tissue) == "Adipose_Subcutaneous",]
stochastic_twas_df_adipose_sub <- stochastic_twas_df[as.character(stochastic_twas_df$tissue) == "Adipose_Subcutaneous",]
fusion_twas_ids <- paste0(as.character(fusion_twas_df_adipose_sub$trait_component),as.character(fusion_twas_df_adipose_sub$gene))
stochastic_twas_ids <- paste0(as.character(stochastic_twas_df_adipose_sub$trait_component),as.character(stochastic_twas_df_adipose_sub$gene))

m <- match(fusion_twas_ids, stochastic_twas_ids)
m_rev <- match(stochastic_twas_ids, fusion_twas_ids)


m.keep = !is.na(m)
fusion_twas_ids = fusion_twas_ids[m.keep]
fusion_twas_df_adipose_sub = fusion_twas_df_adipose_sub[m.keep,]
m_rev.keep = !is.na(m_rev)
stochastic_twas_ids = stochastic_twas_ids[m_rev.keep]
stochastic_twas_df_adipose_sub = stochastic_twas_df_adipose_sub[m_rev.keep,]

#m2 <- match(fusion_twas_ids, stochastic_twas_ids)

#print(sum(as.character(fusion_twas_ids[m2]) != as.character(stochastic_twas_ids)))
#print(sum(as.character(fusion_twas_ids) != as.character(stochastic_twas_ids)))

#print(summary(fusion_twas_df_adipose_sub))

#print(cor(stochastic_twas_df_adipose_sub$nominal_twas_zscore, fusion_twas_df_adipose_sub$nominal_twas_zscore))
#print(cor(stochastic_twas_df_adipose_sub$nominal_stochastic_zscore_mean, fusion_twas_df_adipose_sub$nominal_twas_zscore))
#print(cor(abs(stochastic_twas_df_adipose_sub$nominal_stochastic_zscore_mean -fusion_twas_df_adipose_sub$nominal_twas_zscore), stochastic_twas_df_adipose_sub$nominal_stochastic_zscore_sd))

##############################
# Scatter plot showing correlation of nominal twas PMCES zscores and FUSION twas zscores
##############################
output_file <- paste0(visualize_twas_results_dir, "nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot.pdf")
scatter <- nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot(stochastic_twas_df_adipose_sub, fusion_twas_df_adipose_sub)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")

output_file <- paste0(visualize_twas_results_dir, "nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot_colored_by_stochastic_z_score_sd.pdf")
scatter <- nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot_colored_by_stochastic_z_score_sd(stochastic_twas_df_adipose_sub, fusion_twas_df_adipose_sub)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")


##############################
# Scatter plot showing correlation of absolute difference in nominal twas PMCES zscores and FUSION twas zscores by stochastic z score sd
##############################
output_file <- paste0(visualize_twas_results_dir, "abs_diff_nominal_fusion_zscore_susie_pmces_twas_zscore_by_stochastic_z_score_sd_scatterplot.pdf")
scatter <- abs_diff_nominal_fusion_zscore_susie_pmces_twas_zscore_scatterplot_by_stochastic_z_score_sd(stochastic_twas_df_adipose_sub, fusion_twas_df_adipose_sub)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")



##############################
# Scatter plot showing correlation of nominal twas zscores and mean stochastic twas zscores
##############################
output_file <- paste0(visualize_twas_results_dir, "nominal_zscore_stochastic_mean_zscore_twas_scatterplot.pdf")
scatter <- make_nominal_zscore_stochastic_mean_zscore_scatterplot(stochastic_twas_df)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")

##############################
# Scatter plot showing correlation of nominal twas zscores and mean stochastic twas zscores
##############################
output_file <- paste0(visualize_twas_results_dir, "nominal_zscore_stochastic_mean_zscore_twas_scatterplot_colored_by_vmr.pdf")
scatter <- make_nominal_zscore_stochastic_mean_zscore_scatterplot_colored_by_vmr(stochastic_twas_df)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")


##############################
# Scatter plot showing correlation of nominal twas PIPS and mean stochastic twas PIPS
##############################
output_file <- paste0(visualize_twas_results_dir, "nominal_pip_stochastic_mean_pip_twas_scatterplot.pdf")
scatter <- make_nominal_pip_stochastic_mean_pip_scatterplot(stochastic_twas_df)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")

output_file <- paste0(visualize_twas_results_dir, "nominal_pip_stochastic_mean_pip_twas_scatterplot_colored_by_whole_blood.pdf")
scatter <- make_nominal_pip_stochastic_mean_pip_scatterplot_colored_by_whole_blood(stochastic_twas_df)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")

output_file <- paste0(visualize_twas_results_dir, "nominal_pip_stochastic_mean_pip_twas_scatterplot_colored_by_z_score_sd.pdf")
scatter <- make_nominal_pip_stochastic_mean_pip_scatterplot_colored_by_zscore_sd(stochastic_twas_df)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")

output_file <- paste0(visualize_twas_results_dir, "nominal_pip_stochastic_mean_pip_twas_scatterplot_colored_by_vmr.pdf")
scatter <- make_nominal_pip_stochastic_mean_pip_scatterplot_colored_by_vmr(stochastic_twas_df)
ggsave(scatter, file=output_file, width=7.2, height=5.3, units="in")
