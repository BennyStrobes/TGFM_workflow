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
library(ggrepel)
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

replication_gene_tissue_level_se_barplot_at_single_pip_thresholds <- function(organized_results_dir, pip_threshold) {
	pip_thresholds <- rev(c(pip_threshold))
	
	category_vec <- c()
	count_vec <- c()
	pip_thresh_vec <- c()

	for (pip_iter in 1:length(pip_thresholds)) {

		pip_threshold <- pip_thresholds[pip_iter]

		input_file <- paste0(organized_results_dir, "tgfm_results_hold_out_tissue_summary_", pip_threshold, ".txt")

		raw_df <- read.table(input_file, header=TRUE,sep="\t")


		category <- "no gene-tissue"
		formal_category <- "no gene-tissue"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))


		category <- "same gene, proxy tissue"
		formal_category <- "same gene\n(Spleen)"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

		category <- "same gene, non-proxy tissue"
		formal_category <- "same gene\n(other tissue)"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

		category <- "different gene, tissue"
		formal_category <- "different gene-tissue"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))
	}

	red_color=brewer.pal(n = 9, name = "Reds")[6]
	red_color1=brewer.pal(n = 9, name = "Reds")[3]
	red_color2=brewer.pal(n = 9, name = "Reds")[5]



	ordered_categories <- c("no gene-tissue", "same gene\n(Spleen)", "same gene\n(other tissue)", "different gene-tissue")


	df <- data.frame(counts=count_vec, PIP=factor(pip_thresh_vec,levels=pip_thresholds), gene_category=factor(category_vec, levels=ordered_categories))

	max_count = max(count_vec)
	p <- ggplot(data=df, aes(x=gene_category, y=counts, fill=PIP)) +
  		geom_bar(stat="identity", position="dodge")+
  		figure_theme() +
  		labs(y="\nNo. fine-mapped (PIP > 0.5)\ngene-tissue pairs", x="") +
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  		scale_fill_manual(values=rev(c(red_color)))+
  		geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25) +
  		ylim(0.0,max_count+5) +
  		theme(plot.title = element_text(hjust = 0.5)) + 
  		theme(legend.position="none")

  	return(p)


}


replication_gene_tissue_level_se_barplot_at_various_pip_thresholds <- function(organized_results_dir, version="replication") {
	pip_thresholds <- rev(c("0.3", "0.4", "0.5"))
	
	category_vec <- c()
	count_vec <- c()
	pip_thresh_vec <- c()

	for (pip_iter in 1:length(pip_thresholds)) {

		pip_threshold <- pip_thresholds[pip_iter]

		input_file <- paste0(organized_results_dir, "tgfm_results_hold_out_tissue_summary_", pip_threshold, ".txt")
		if (version=="subsampling") {
			input_file <- paste0(organized_results_dir, "tgfm_results_whole_blood_subsampled_tissue_summary_", pip_threshold, ".txt")
		}
		raw_df <- read.table(input_file, header=TRUE,sep="\t")


		category <- "no gene-tissue"
		formal_category <- "no gene-tissue"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

		category <- "same gene, replication tissue"
		if (version == "replication") {
			formal_category <- "same gene\n(PBMC)"
		}
		if (version == "subsampling") {
			formal_category <- "same gene\n(Whole Blood subsampled)"
		}
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

		category <- "same gene, proxy tissue"
		formal_category <- "same gene\n(Spleen)"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

		category <- "same gene, non-proxy tissue"
		formal_category <- "same gene\n(other tissue)"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

		category <- "different gene, tissue"
		formal_category <- "different gene-tissue"
		counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
		category_vec <- c(category_vec, formal_category)
		count_vec <- c(count_vec, counts)
		pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))
	}

	red_color=brewer.pal(n = 9, name = "Reds")[7]
	red_color1=brewer.pal(n = 9, name = "Reds")[3]
	red_color2=brewer.pal(n = 9, name = "Reds")[5]


	if (version == "replication") {
		ordered_categories <- c("no gene-tissue", "same gene\n(PBMC)", "same gene\n(Spleen)", "same gene\n(other tissue)", "different gene-tissue")
	}
	if (version == "subsampling") {
		ordered_categories <- c("no gene-tissue", "same gene\n(Whole Blood subsampled)", "same gene\n(Spleen)", "same gene\n(other tissue)", "different gene-tissue")

	}

	df <- data.frame(counts=count_vec, PIP=factor(pip_thresh_vec,levels=pip_thresholds), gene_category=factor(category_vec, levels=ordered_categories))

	max_count = max(count_vec)
	p <- ggplot(data=df, aes(x=gene_category, y=counts, fill=PIP)) +
  		geom_bar(stat="identity", position="dodge")+
  		figure_theme() +
  		labs(y="No. fine-mapped\ngene-tissue pairs", x="", title=paste0("Whole Blood ", version)) +
  		theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  		scale_fill_manual(values=rev(c(red_color1, red_color2, red_color)))+
  		geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25) +
  		ylim(0.0,max_count+5) +
  		theme(plot.title = element_text(hjust = 0.5))

  	return(p)


}

get_gene_tissue_level_se_barplot_at_single_pip_threshold_data <- function(tgfm_organized_results_dir, pip_threshold) {
	pip_thresholds <- rev(c(pip_threshold))
	
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
	formal_category <- "same gene (proxy tissue)"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, formal_category)
	count_vec <- c(count_vec, counts)
	pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

	category <- "same gene, non-proxy tissue"
	formal_category <- "same gene (non-proxy tissue)"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, formal_category)
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


	df <- data.frame(counts=count_vec, PIP=factor(pip_thresh_vec,levels=pip_thresholds), gene_category=factor(category_vec, levels=c("no gene-tissue", "same gene (proxy tissue)", "same gene (non-proxy tissue)", "different gene, tissue")))

	return(df)
}



gene_tissue_level_se_barplot_at_single_pip_threshold <- function(tgfm_organized_results_dir, pip_threshold) {
	pip_thresholds <- rev(c(pip_threshold))
	
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
	formal_category <- "same gene\n(proxy tissue)"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, formal_category)
	count_vec <- c(count_vec, counts)
	pip_thresh_vec <- c(pip_thresh_vec, rep(pip_threshold, length(counts)))

	category <- "same gene, non-proxy tissue"
	formal_category <- "same gene\n(non-proxy tissue)"
	counts = sum(as.character(raw_df$gene_tissue_level_description) == category)
	category_vec <- c(category_vec, formal_category)
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


	df <- data.frame(counts=count_vec, PIP=factor(pip_thresh_vec,levels=pip_thresholds), gene_category=factor(category_vec, levels=c("no gene-tissue", "same gene\n(proxy tissue)", "same gene\n(non-proxy tissue)", "different gene, tissue")))

	p <- ggplot(data=df, aes(x=gene_category, y=counts, fill=PIP)) +
  		geom_bar(stat="identity", position="dodge")+
  		figure_theme() +
  		labs(y="\nNo. fine-mapped\ngene-tissue pairs", x="") +
  		theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  		scale_fill_manual(values=rev(c(red_color)))+
  		geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25) +
  		ylim(0.0,105.0) + 
  		theme(legend.position="none")

  	return(p)
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

make_ablated_hit_summary_delta_pip_stack_barplot_cross_trait <- function(df_raw) {
	new_vec <- paste0(df_raw$trait_name, ":", df_raw$original_gt_pair)

	df_raw$delta_nm_var_pip[df_raw$delta_nm_var_pip <= 0] = 0.0
	df_raw$delta_gene_pip[df_raw$delta_gene_pip <= 0] = 0.0
	
	delta_pip_vec <- c(df_raw$delta_nm_var_pip, df_raw$delta_gene_pip)
	gt_pair_vec <- c(new_vec, new_vec)
	pip_type_vec <- c(rep("NM-Variant", length(df_raw$delta_nm_var_pip)), rep("Gene-Tissue", length(df_raw$delta_gene_pip)))

	df <- data.frame(delta_pip=delta_pip_vec, gt_pair=gt_pair_vec, pip_type=pip_type_vec)

	totes = df_raw$delta_nm_var_pip + df_raw$delta_gene_pip

	indexes <- order(totes)

	df$gt_pair = factor(df$gt_pair,levels=as.character(new_vec)[indexes])


	pp <- ggplot(df, aes(fill=pip_type_vec, y=delta_pip, x=gt_pair)) + 
    	geom_bar(position="stack", stat="identity") +
    	figure_theme() +
    	theme(axis.text.x=element_blank()) +
    	labs(x="",y="Delta PIP", fill="") +
    	scale_fill_manual(values=c("#999999", "#56B4E9")) + 
    	theme(legend.position="bottom")
    return(pp)


}


make_ablated_hit_summary_delta_pip_stack_barplot_cross_trait_v2 <- function(df_v2) {
	## V2 proxy gene-tissue
	new_vec2 <- paste0(df_v2$trait_name, ":", df_v2$original_gt_pair)
	df_v2$delta_nm_var_pip[df_v2$delta_nm_var_pip <= 0] = 0.0
	df_v2$delta_nontagging_gene_pip[df_v2$delta_nontagging_gene_pip <= 0] = 0.0
	df_v2$delta_tagging_gene_pip[df_v2$delta_tagging_gene_pip <= 0] = 0.0
	
	delta_pip_vec2 <- c(df_v2$delta_nm_var_pip, df_v2$delta_nontagging_gene_pip, df_v2$delta_tagging_gene_pip)
	gt_pair_vec2 <- c(new_vec2, new_vec2, new_vec2)
	pip_type_vec2 <- c(rep("NM-Variant", length(df_v2$delta_nm_var_pip)), rep("Non-Proxy Gene-Tissue", length(df_v2$delta_nontagging_gene_pip)), rep("Proxy Gene-Tissue", length(df_v2$delta_tagging_gene_pip)))

	df2 <- data.frame(delta_pip=delta_pip_vec2, gt_pair=gt_pair_vec2, pip_type=pip_type_vec2)

	totes2 = df_v2$delta_nm_var_pip + df_v2$delta_nontagging_gene_pip + df_v2$delta_tagging_gene_pip
	indexes2 <- order(totes2)
	df2$gt_pair = factor(df2$gt_pair,levels=as.character(new_vec2)[indexes2])
	ppp2 <- ggplot(df2, aes(fill=pip_type, y=delta_pip, x=gt_pair)) + 
    	geom_bar(position="stack", stat="identity") +
    	figure_theme() +
    	theme(axis.text.x=element_blank()) +
    	labs(x="",y="Delta PIP", fill="") +
    	scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
    	theme(legend.position="bottom")

}


make_ablated_hit_summary_delta_pip_stack_barplot_cross_trait_tagging_tissue_strat <- function(df_raw) {

	## V1 no proxy gene-tissue
	df_v1 = df_raw[df_raw$delta_tagging_gene_pip==0,]

	new_vec <- paste0(df_v1$trait_name, ":", df_v1$original_gt_pair)
	df_v1$delta_nm_var_pip[df_v1$delta_nm_var_pip <= 0] = 0.0
	df_v1$delta_nontagging_gene_pip[df_v1$delta_nontagging_gene_pip <= 0] = 0.0
	
	delta_pip_vec <- c(df_v1$delta_nm_var_pip, df_v1$delta_nontagging_gene_pip)
	gt_pair_vec <- c(new_vec, new_vec)
	pip_type_vec <- c(rep("NM-Variant", length(df_v1$delta_nm_var_pip)), rep("Non-Proxy Gene-Tissue", length(df_v1$delta_nontagging_gene_pip)))

	df <- data.frame(delta_pip=delta_pip_vec, gt_pair=gt_pair_vec, pip_type=pip_type_vec)

	totes = df_v1$delta_nm_var_pip + df_v1$delta_nontagging_gene_pip
	indexes <- order(totes)
	df$gt_pair = factor(df$gt_pair,levels=as.character(new_vec)[indexes])
	ppp1 <- ggplot(df, aes(fill=pip_type, y=delta_pip, x=gt_pair)) + 
    	geom_bar(position="stack", stat="identity") +
    	figure_theme() +
    	theme(axis.text.x=element_blank()) +
    	labs(x="",y="Delta PIP", fill="", title="Ablated gene-tissue hax no proxy gene-tissue") +
    	scale_fill_manual(values=c("#999999", "#56B4E9")) +
    	theme(legend.position="bottom")




	## V2 proxy gene-tissue
	df_v2 = df_raw[df_raw$delta_tagging_gene_pip!=0,]

	new_vec2 <- paste0(df_v2$trait_name, ":", df_v2$original_gt_pair)
	df_v2$delta_nm_var_pip[df_v2$delta_nm_var_pip <= 0] = 0.0
	df_v2$delta_nontagging_gene_pip[df_v2$delta_nontagging_gene_pip <= 0] = 0.0
	df_v2$delta_tagging_gene_pip[df_v2$delta_tagging_gene_pip <= 0] = 0.0
	
	delta_pip_vec2 <- c(df_v2$delta_nm_var_pip, df_v2$delta_nontagging_gene_pip, df_v2$delta_tagging_gene_pip)
	gt_pair_vec2 <- c(new_vec2, new_vec2, new_vec2)
	pip_type_vec2 <- c(rep("NM-Variant", length(df_v2$delta_nm_var_pip)), rep("Non-Proxy Gene-Tissue", length(df_v2$delta_nontagging_gene_pip)), rep("Proxy Gene-Tissue", length(df_v2$delta_tagging_gene_pip)))

	df2 <- data.frame(delta_pip=delta_pip_vec2, gt_pair=gt_pair_vec2, pip_type=pip_type_vec2)

	totes2 = df_v2$delta_nm_var_pip + df_v2$delta_nontagging_gene_pip + df_v2$delta_tagging_gene_pip
	indexes2 <- order(totes2)
	df2$gt_pair = factor(df2$gt_pair,levels=as.character(new_vec2)[indexes2])
	ppp2 <- ggplot(df2, aes(fill=pip_type, y=delta_pip, x=gt_pair)) + 
    	geom_bar(position="stack", stat="identity") +
    	figure_theme() +
    	theme(axis.text.x=element_blank()) +
    	labs(x="",y="Delta PIP", fill="", title="Ablated gene-tissue has proxy gene-tissue") +
    	scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
    	theme(legend.position="bottom")

    legender <- get_legend(ppp2)
    joint <- plot_grid(ppp1+ theme(legend.position="None"),ppp2 + theme(legend.position="None"),ncol=2)

    pp <- plot_grid(joint,legender,ncol=1, rel_heights=c(1,.15))

    return(pp)


}

get_replication_analysis_histogram_data <-  function(repication_data_file, replication_name, original_tissue, replicating_tissue, replicating_tissue_name_readable, titler) {
	df <- read.table(repication_data_file, header=TRUE, sep="\t")

	all_tissue_names <- as.character(unique(df$new_tissue))

	# Remore replicating tissue
	other_tissue_names <- all_tissue_names[all_tissue_names!=replicating_tissue]


	tissue_names_arr <- c()
	avg_pip_arr <- c()


	# loop Through tissues
	for (tissue_iter in 1:length(all_tissue_names)) {
		tissue_name <- all_tissue_names[tissue_iter]
		avg_pip = mean(df$new_tissue_pip[as.character(df$new_tissue) == tissue_name])

		tissue_names_arr <- c(tissue_names_arr, tissue_name)
		avg_pip_arr <- c(avg_pip_arr, avg_pip)

	}

	df2 <- data.frame(PIP=avg_pip_arr, tissue=tissue_names_arr)

	return(df2)


}

make_replication_analysis_histogram2 <- function(repication_data_file, replication_name, original_tissue, replicating_tissue, replicating_tissue_name_readable, titler) {
	df <- read.table(repication_data_file, header=TRUE, sep="\t")

	all_tissue_names <- as.character(unique(df$new_tissue))

	# Remore replicating tissue
	other_tissue_names <- all_tissue_names[all_tissue_names!=replicating_tissue]


	tissue_names_arr <- c()
	avg_pip_arr <- c()


	# loop Through tissues
	for (tissue_iter in 1:length(other_tissue_names)) {
		tissue_name <- other_tissue_names[tissue_iter]
		avg_pip = mean(df$new_tissue_pip[as.character(df$new_tissue) == tissue_name])

		tissue_names_arr <- c(tissue_names_arr, tissue_name)
		avg_pip_arr <- c(avg_pip_arr, avg_pip)

	}

	df2 <- data.frame(PIP=avg_pip_arr, tissue=tissue_names_arr)


	replication_avg_pip = mean(df$new_tissue_pip[as.character(df$new_tissue) == replicating_tissue])

	df3 <- data.frame(PIPPY=c(replication_avg_pip), labeler=c(replicating_tissue_name_readable), density=c(10))

	red_color=brewer.pal(n = 9, name = "Reds")[7]


	p <- ggplot(df2, aes(PIP)) +
		geom_histogram(fill="skyblue",color="grey15", alpha=0.9, binwidth=0.02, boundary = 0) + 
  		figure_theme() +
  		geom_vline(xintercept = replication_avg_pip,color = red_color, size=1.5) + 
  		#xlim(-.05, .35) +
		geom_text_repel(data=df3, aes(x=PIPPY, y=density, label=labeler), color=red_color,size=4.4, nudge_x=-.02) +
		labs(y="\n\nNo. tissues", x="Average gene-tissue PIP / tissue") +
		theme(plot.title = element_text(hjust = 0.5))

  	return(p)
}

make_replication_analysis_histogram <- function(repication_data_file, replication_name, original_tissue, replicating_tissue, replicating_tissue_name_readable, titler) {
	df <- read.table(repication_data_file, header=TRUE, sep="\t")

	all_tissue_names <- as.character(unique(df$new_tissue))

	# Remore replicating tissue
	other_tissue_names <- all_tissue_names[all_tissue_names!=replicating_tissue]


	tissue_names_arr <- c()
	avg_pip_arr <- c()


	# loop Through tissues
	for (tissue_iter in 1:length(other_tissue_names)) {
		tissue_name <- other_tissue_names[tissue_iter]
		avg_pip = mean(df$new_tissue_pip[as.character(df$new_tissue) == tissue_name])

		tissue_names_arr <- c(tissue_names_arr, tissue_name)
		avg_pip_arr <- c(avg_pip_arr, avg_pip)

	}

	df2 <- data.frame(PIP=avg_pip_arr, tissue=tissue_names_arr)


	replication_avg_pip = mean(df$new_tissue_pip[as.character(df$new_tissue) == replicating_tissue])

	df3 <- data.frame(PIPPY=c(replication_avg_pip), labeler=c(replicating_tissue_name_readable), density=c(10))

		red_color=brewer.pal(n = 9, name = "Reds")[7]


	p <- ggplot(df2, aes(PIP)) +
  		geom_density(fill="skyblue", alpha=0.6) +
  		figure_theme() +
  		geom_vline(xintercept = replication_avg_pip,color = red_color, size=1.5) + 
  		xlim(0, .35) +
		geom_text_repel(data=df3, aes(x=PIPPY, y=density, label=labeler), color=red_color,size=4.4, nudge_x=-.02) +
		labs(y="\n\nNo. tissues", x="Average gene-tissue PIP / tissue") +
		theme(plot.title = element_text(hjust = 0.5))

  	return(p)
}


tgfm_organized_results_dir <- args[1]
sc_pb_tgfm_organized_results_dir <- args[2]
trait_names_file <- args[3]
trait_names_file2 <- args[4]
output_dir <- args[5]
tissue_replication_results_dir <- args[6]


print(tissue_replication_results_dir)


trait_df <- read.table(trait_names_file, header=TRUE,sep="\t")
trait_names <- as.character(trait_df$study_name)
hold_out_tissues <- as.character(trait_df$tissue_remove)



trait_df2 <- read.table(trait_names_file2, header=TRUE,sep="\t")
trait_names_wb <- as.character(trait_df2$study_name)



ablated_hit_summary_file = paste0(tgfm_organized_results_dir, "tgfm_results_hold_out_tissue_ablated_hit_summary.txt")
ablated_hit_df <- read.table(ablated_hit_summary_file,header=TRUE,sep="\t")





##################################################
# Make replication analysis histogram
##################################################
# PBMC replication
replication_name="Whole_Blood_PBMC"
original_tissue="Whole_Blood"
replicating_tissue="PBMC"
repication_data_file <- paste0(tissue_replication_results_dir, replication_name, "_replication_pip_0.5_raw_replication_results.txt")
rep_histo <- make_replication_analysis_histogram(repication_data_file, replication_name, original_tissue, replicating_tissue , "PBMC", "Whole Blood replication")
output_file = paste0(output_dir, replication_name, "_replication_density.pdf")
ggsave(rep_histo, file=output_file, width=7.2, height=4.1, units="in")

# PBMC replication
replication_name="Whole_Blood_PBMC"
original_tissue="Whole_Blood"
replicating_tissue="PBMC"
repication_data_file <- paste0(tissue_replication_results_dir, replication_name, "_replication_pip_0.5_raw_replication_results.txt")
rep_histo2 <- make_replication_analysis_histogram2(repication_data_file, replication_name, original_tissue, replicating_tissue , "PBMC", "Whole Blood replication")
output_file = paste0(output_dir, replication_name, "_replication_histogram.pdf")
ggsave(rep_histo2, file=output_file, width=7.2, height=4.1, units="in")


# Whole blood subsampled replication
replication_name="Whole_Blood_Whole_Blood_subsampled"
original_tissue="Whole_Blood"
replicating_tissue="Whole_Blood_subsampled"
repication_data_file <- paste0(tissue_replication_results_dir, replication_name, "_replication_pip_0.5_raw_replication_results.txt")
subsampled_histo <- make_replication_analysis_histogram(repication_data_file, replication_name, original_tissue, replicating_tissue, "Whole Blood\n(subsampled)", "Whole Blood subsampling")
output_file = paste0(output_dir, replication_name, "_replication_density.pdf")
ggsave(subsampled_histo, file=output_file, width=7.2, height=4.1, units="in")


# Whole blood subsampled replication
replication_name="Whole_Blood_Whole_Blood_subsampled"
original_tissue="Whole_Blood"
replicating_tissue="Whole_Blood_subsampled"
repication_data_file <- paste0(tissue_replication_results_dir, replication_name, "_replication_pip_0.5_raw_replication_results.txt")
subsampled_histo <- make_replication_analysis_histogram2(repication_data_file, replication_name, original_tissue, replicating_tissue, "Whole Blood\n(subsampled)", "Whole Blood subsampling")
output_file = paste0(output_dir, replication_name, "_replication_histogram.pdf")
ggsave(subsampled_histo, file=output_file, width=7.2, height=4.1, units="in")

#######################
# Make se barplot showing number of gene-tissue pairs discovered in the ablated analyisis (for loci where PIP > 0.5 in the original analysis)
# at single pip threshold
#######################
pip_threshold = 0.5
n_ablation_hits_bar <- gene_tissue_level_se_barplot_at_single_pip_threshold(tgfm_organized_results_dir, pip_threshold) + theme(plot.margin = margin(5.5, 5.5, 0.0, 5.5, "points"))
output_file <- paste0(output_dir, "held_out_tissue_gene_tissue_level_summary_barplot_pip_", pip_threshold, ".pdf")
ggsave(n_ablation_hits_bar, file=output_file, width=7.2, height=4.6, units="in")

###########################
# Make figure 5 with cowplot
##########################
#fig_5bc <- plot_grid(subsampled_histo, rep_histo, ncol=2, labels=c("b", "c"))
#rep_histo <- rep_histo + labs(x="Average replicating gene−tissue PIP / tissue", y="\nNo. tissues")
n_ablation_hits_bar = n_ablation_hits_bar + theme(plot.margin = margin(5.5, 5.5, 0.0, 5.5, "points")) + labs(x=NULL)
rep_histo <- plot_grid(NULL, rep_histo2, ncol=2, rel_widths=c(.012, 1.0))
fig_5 <- plot_grid(n_ablation_hits_bar, rep_histo, ncol=1, rel_heights=c(1.5,1), labels=c("a","b"))
output_file <- paste0(output_dir, "figure5.pdf")
ggsave(fig_5, file=output_file, width=7.2, height=4.0, units="in")


###########################
# Save figure 5 data to supplementary table
##########################
# 5a supptable
if (FALSE) {
supp_table_file = paste0(output_dir, "suppTable_figure5a_numerical.txt")
pip_threshold = 0.5
supp_table_df <- get_gene_tissue_level_se_barplot_at_single_pip_threshold_data(tgfm_organized_results_dir, pip_threshold)
write.table(supp_table_df, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
}

# 5b supptable
if (FALSE) {
supp_table_file = paste0(output_dir, "suppTable_figure5b_numerical.txt")
pip_threshold = 0.5
replication_name="Whole_Blood_PBMC"
original_tissue="Whole_Blood"
replicating_tissue="PBMC"
repication_data_file <- paste0(tissue_replication_results_dir, replication_name, "_replication_pip_0.5_raw_replication_results.txt")
supp_table_df <- get_replication_analysis_histogram_data(repication_data_file, replication_name, original_tissue, replicating_tissue , "PBMC", "Whole Blood replication")
write.table(supp_table_df, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
print(supp_table_file)
}

#######################
# Make se barplot showing number of gene-tissue pairs discovered in the replication analyisis (for loci where PIP > 0.5 in the original analysis)
# for various pip thresholds
#######################
if (FALSE) {
pp <- replication_gene_tissue_level_se_barplot_at_various_pip_thresholds(sc_pb_tgfm_organized_results_dir, version="replication")
output_file <- paste0(output_dir, "replication_tissue_gene_tissue_level_summary_barplot_various_pip.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")

print(sc_pb_tgfm_organized_results_dir)
pp <- replication_gene_tissue_level_se_barplot_at_various_pip_thresholds(tgfm_organized_results_dir, version="subsampling")
output_file <- paste0(output_dir, "subsampling_tissue_gene_tissue_level_summary_barplot_various_pip.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")
}
#####################
#Make stacked bar plot showing delta nm variant and delta gene after ablation
####################
if (FALSE) {
pp <- make_ablated_hit_summary_delta_pip_stack_barplot_cross_trait(ablated_hit_df)
output_file <- paste0(output_dir, "ablated_hit_delta_pip_stacked_barplot_cross_trait.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")

pp <- make_ablated_hit_summary_delta_pip_stack_barplot_cross_trait_v2(ablated_hit_df)
output_file <- paste0(output_dir, "ablated_hit_delta_pip_stacked_barplot_cross_trait_v2.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")

pp <- make_ablated_hit_summary_delta_pip_stack_barplot_cross_trait_tagging_tissue_strat(ablated_hit_df)
output_file <- paste0(output_dir, "ablated_hit_delta_pip_stacked_barplot_cross_trait_tagging_tiss_strat.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")
}

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
if (FALSE) {
pp <- gene_tissue_level_se_barplot_at_various_pip_thresholds(tgfm_organized_results_dir)
output_file <- paste0(output_dir, "held_out_tissue_gene_tissue_level_summary_barplot_various_pip.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.6, units="in")
}










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


