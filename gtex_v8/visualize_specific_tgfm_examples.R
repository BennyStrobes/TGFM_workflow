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
library(ggnewscale)
options(warn=1)
options(bitmapType='cairo')


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



make_tgfm_manhatten_plot_for_given_window <- function(snp_df, gene_df, link_df, trait_name, ensamble_tissue_name, gene_tissue_name, window_name) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	# Split both snps and genes into hits and not
	snp_hits_df = snp_df[snp_df$TGFM_PIP >= .2,]
	snp_null_df = snp_df[snp_df$TGFM_PIP < .2,]
	gene_hits_df = gene_df[gene_df$TGFM_PIP >= .2,]
	gene_null_df = gene_df[gene_df$TGFM_PIP < .2,]



	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]

	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.2) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.2) + 
		  scale_colour_gradient(low = "thistle1", high = "red") +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.9) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=1.9) + 
		  #geom_text_repel(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=gene_id_tissue_name, color=TGFM_PIP),size=2.8) +
		  #geom_errorbar(data=gene_hits_df, aes(x=gene_tss_mb,y=gwas_neg_log10_p_mean, ymax=gwas_neg_log10_p_ub, ymin=gwas_neg_log10_p_lb,color=TGFM_PIP),size=.3, width=.1) +
		  scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue3") +
		  figure_theme() +
		  labs(x=paste0("Position [MB]\nChromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name) + 
		  theme(legend.position="bottom")

	return(pp)
}




















##################
# Load in data
##################
specific_examples_input_file = args[1]
visualize_specific_tgfm_examples_dir = args[2]

# Load in specific examples data
specific_examples_data <- read.table(specific_examples_input_file, header=TRUE, sep="\t")



# Loop through specific examples
n_examples = dim(specific_examples_data)[1]
for (example_iter in 1:n_examples) {
	# Extract relevent info for this example
	trait_name <- specific_examples_data$trait[example_iter]
	ensamble_tissue_name <- specific_examples_data$gene_tissue[example_iter]
	gene_id <- specific_examples_data$gene_id[example_iter]
	tissue_name <- specific_examples_data$tissue[example_iter]
	window_name <- specific_examples_data$window_name[example_iter]
	gene_tissue_name <- paste0(gene_id, "-", tissue_name)
	# Input files
	snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
	gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
	link_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_link_df.txt")
	# Load in input
	snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
	gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")
	#link_df <- read.table(link_df_input_file, header=TRUE, sep="\t")
	link_df=0

	# Make manhatten plot for a given window/trait
	tgfm_manhatten_plot <- make_tgfm_manhatten_plot_for_given_window(snp_df, gene_df, link_df, trait_name, ensamble_tissue_name, gene_tissue_name, window_name)
	output_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_tgfm_manhattan.pdf")
	ggsave(tgfm_manhatten_plot, file=output_file, width=7.2, height=4.0, units="in")

}








