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



gene_level_se_barplot <- function(raw_df) {
	category_vec <- c()
	count_vec <- c()

	category <- "no gene"
	counts = sum(as.character(raw_df$gene_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)

	category <- "same gene"
	counts = sum(as.character(raw_df$gene_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)

	category <- "different gene"
	counts = sum(as.character(raw_df$gene_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)


	df <- data.frame(counts=count_vec, gene_category=factor(category_vec, levels=c("no gene", "same gene", "different gene")))

	p <- ggplot(data=df, aes(x=gene_category, y=counts)) +
  		geom_bar(stat="identity", fill="steelblue")+
  		figure_theme() +
  		labs(y="No. fine-mapped (PIP > 0.5) genes", x="") +
  		geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25)

  	return(p)
}

gene_tissue_level_se_barplot <- function(raw_df) {
	category_vec <- c()
	count_vec <- c()

	category <- "no gene-tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)

	category <- "same gene, tagging tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)

	category <- "same gene, non-tagging tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)

	category <- "different gene, tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)



	df <- data.frame(counts=count_vec, gene_category=factor(category_vec, levels=c("no gene-tissue", "same gene, tagging tissue", "same gene, non-tagging tissue", "different gene, tissue")))

	p <- ggplot(data=df, aes(x=gene_category, y=counts)) +
  		geom_bar(stat="identity", fill="steelblue")+
  		figure_theme() +
  		labs(y="No. fine-mapped (PIP > 0.5)\ngene-tissue pairs", x="") +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25) +
  		ylim(0.0,105.0)
  	joint <- plot_grid(NULL, p,ncol=1,rel_heights=c(.06,1))

  	return(p)
}

make_before_after_scatterplot <- function(summary_file, trait_name, held_out_tissue) {
	df <- read.table(summary_file, header=TRUE, sep="\t")

	pp <- ggplot(df, aes(x=gene_tissue_pip, y=original_pip,color=tagging_tissue)) + geom_point() +
		figure_theme() + 
		labs(x="Gene-Tissue PIP (tissue subset analysis)", y="Gene-Tissue PIP (full analysis)", title=paste0(trait_name, " / ", held_out_tissue)) +
		geom_abline(slope=1)
	return(pp)
}

make_before_after_scatterplot_cross_trait <- function(tgfm_organized_results_dir, trait_names, hold_out_tissues) {
	gt_pip_vec <- c()
	gt_pip_orig_vec <- c()
	tagging_vec <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		held_out_tissue <- hold_out_tissues[trait_iter]
		summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_ignore_", held_out_tissue, "_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_gene_tissue_hit_summary.txt")

		df <- read.table(summary_file, header=TRUE, sep="\t")

		gt_pip_vec <- c(gt_pip_vec, df$gene_tissue_pip)
		gt_pip_orig_vec <- c(gt_pip_orig_vec, df$original_pip)
		tagging_vec <- c(tagging_vec, as.character(df$tagging_tissue))

	}

	new_df <- data.frame(gene_tissue_pip=gt_pip_vec, original_pip=gt_pip_orig_vec, tagging_tissue=factor(tagging_vec))

	print(dim(new_df))

	pp <- ggplot(new_df, aes(x=gene_tissue_pip, y=original_pip,color=tagging_tissue)) + geom_point() +
		figure_theme() + 
		labs(x="Gene-Tissue PIP (tissue subset analysis)", y="Gene-Tissue PIP (full analysis)",color="") +
		geom_abline(slope=1)
	return(pp)


}

make_ablated_hit_summary_scatterplot_cross_trait <- function(df, title="") {
	df$gt_pair_from_proxy_tissue = factor(df$gt_pair_from_proxy_tissue, levels=c("proxy tissue", "other tissue"))

	pp <- ggplot(df, aes(y=ablated_gt_pair_pip, x=original_gt_pair_pip,color=gt_pair_from_proxy_tissue)) + geom_point(size=1.7) +
		figure_theme() + 
		labs(y="Gene-Tissue PIP (tissue ablated analysis)", x="Gene-Tissue PIP (full analysis)",color="",title=title) +
		geom_abline(slope=1) +
		ylim(0,1) + 
		theme(legend.position="bottom")

	if (length(unique(df$gt_pair_from_proxy_tissue)) == 2) {
		pp = pp + scale_color_manual(values=c('#56B4E9', '#999999'))
	} else {
		pp = pp + scale_color_manual(values=c('#999999'))
	}
	return(pp)
}

gene_tissue_level_se_barplot_at_various_pip_thresholds <- function(tgfm_organized_results_dir) {
	pip_thresholds <- rev(c("0.3", "0.4", "0.5"))
	
	category_vec <- c()
	count_vec <- c()
	pip_thresh_vec <- c()

	for (pip_iter in 1:length(pip_thresholds)) {

	pip_threshold <- pip_thresholds[pip_iter]

	input_file <- paste0(tgfm_organized_results_dir, "tgfm_results_hold_out_tissue_summary_", pip_threshold, ".txt")
	raw_df <- read.table(input_file, header=TRUE,sep="\t")

	category <- "no gene-tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)
	pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

	category <- "same gene, proxy tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)
	pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

	category <- "same gene, non-proxy tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)
	pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

	category <- "different gene, tissue"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, category)
	count_vec <- c(count_vec, counts)
	pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

	}

	red_color=brewer.pal(n = 9, name = "Reds")[7]
	red_color1=brewer.pal(n = 9, name = "Reds")[3]
	red_color2=brewer.pal(n = 9, name = "Reds")[5]


	df <- data.frame(counts=count_vec, PIP=factor(pip_thresh_vec,levels=pip_thresholds), gene_category=factor(category_vec, levels=c("no gene-tissue", "same gene, proxy tissue", "same gene, non-proxy tissue", "different gene, tissue")))

	p <- ggplot(data=df, aes(x=gene_category, y=counts, fill=PIP)) +
  		geom_bar(stat="identity", position="dodge")+
  		figure_theme() +
  		labs(y="No. fine-mapped (PIP > 0.5)\ngene-tissue pairs", x="") +
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  		scale_fill_manual(values=rev(c(red_color1, red_color2, red_color)))+
  		geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25) +
  		ylim(0.0,105.0)

  	return(p)
}



tgfm_organized_results_dir <- args[1]
trait_names_file <- args[2]
output_dir <- args[3]


trait_df <- read.table(trait_names_file, header=TRUE,sep="\t")
trait_names <- as.character(trait_df$study_name)
hold_out_tissues <- as.character(trait_df$tissue_remove)


if (FALSE) {
#####################
# Make scatter plot showing change in gene-tissue pair pips between ablated analysis and non ablated analysis (for loci where PIP > 0.5 in original analysis)
####################
ablated_hit_summary_file = paste0(tgfm_organized_results_dir, "tgfm_results_hold_out_tissue_ablated_hit_summary.txt")
ablated_hit_df <- read.table(ablated_hit_summary_file,header=TRUE,sep="\t")
pp <- make_ablated_hit_summary_scatterplot_cross_trait(ablated_hit_df)
output_file <- paste0(output_dir, "ablated_hit_summary_gt_pip_scatter.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")


#####################
# Make scatter plot showing change in gene-tissue pair pips between ablated analysis and non ablated analysis (for loci where PIP > 0.5 in original analysis)
# But trait stratefied
####################
ablated_hit_summary_file = paste0(tgfm_organized_results_dir, "tgfm_results_hold_out_tissue_ablated_hit_summary.txt")
ablated_hit_df <- read.table(ablated_hit_summary_file,header=TRUE,sep="\t")
pp1 <- make_ablated_hit_summary_scatterplot_cross_trait(ablated_hit_df[as.character(ablated_hit_df$trait_has_proxy_tissue) == "trait_has_proxy_tissue",], title="Traits where\nablated tissue has proxy tissue")
legender = get_legend(pp1)
pp2 <- make_ablated_hit_summary_scatterplot_cross_trait(ablated_hit_df[as.character(ablated_hit_df$trait_has_proxy_tissue) == "trait_has_no_proxy_tissue",], title="Traits where\nablated tissue has no proxy tissue")
joint_pp <- plot_grid(pp2 + theme(legend.position="none"),pp1 + theme(legend.position="none"),ncol=2)
joint_pp2 <- plot_grid(joint_pp, legender, ncol=1, rel_heights=c(1,.067))
output_file <- paste0(output_dir, "ablated_hit_summary_gt_pip_scatter_trait_stratefied.pdf")
ggsave(joint_pp2, file=output_file, width=7.2, height=4.1, units="in")
}



#######################
# Make se barplot showing number of gene-tissue pairs discovered in the ablated analyisis (for loci where PIP > 0.5 in the original analysis)
# for various pip thresholds
#######################
pp <- gene_tissue_level_se_barplot_at_various_pip_thresholds(tgfm_organized_results_dir)
output_file <- paste0(output_dir, "held_out_tissue_gene_tissue_level_summary_barplot_various_pip.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")











if (FALSE) {
for (trait_iter in 1:length(trait_names)) {
	trait_name <- trait_names[trait_iter]
	held_out_tissue <- hold_out_tissues[trait_iter]
	summary_file <- paste0(tgfm_organized_results_dir, "tgfm_results_", trait_name, "_component_gene_ignore_", held_out_tissue, "_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_gene_tissue_hit_summary.txt")
	pp <- make_before_after_scatterplot(summary_file, trait_name, held_out_tissue)
	output_file <- paste0(output_dir, "before_after_scatter_", trait_name, ".pdf")
	ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")
}
}

if (FALSE) {
pp <- make_before_after_scatterplot_cross_trait(tgfm_organized_results_dir, trait_names, hold_out_tissues)
output_file <- paste0(output_dir, "before_after_scatter_cross_trait.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")
}


if (FALSE) {
held_out_tissue_df <- read.table(held_out_tissue_summary_file, header=TRUE,sep="\t")
}
if (FALSE) {
#################
# Make se barplot showing gene level description
pp <- gene_level_se_barplot(held_out_tissue_df)
output_file <- paste0(output_dir, "held_out_tissue_gene_level_summary_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")

#################
# Make se barplot showing gene-tissue level description
pp <- gene_tissue_level_se_barplot(held_out_tissue_df)
output_file <- paste0(output_dir, "held_out_tissue_gene_tissue_level_summary_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")
}


