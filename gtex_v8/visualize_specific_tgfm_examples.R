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
library(RColorBrewer)
options(warn=1)
options(bitmapType='cairo')


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


poster_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=19), text = element_text(size=19),axis.text=element_text(size=19), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=13), legend.title = element_text(size=18)))
}



make_tgfm_manhatten_plot_for_given_window <- function(snp_df, gene_df, link_df, trait_name, ensamble_tissue_name, gene_tissue_name, window_name) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	# Split both snps and genes into hits and not
	snp_hits_df = snp_df[snp_df$TGFM_PIP >= .5,]
	snp_null_df = snp_df[snp_df$TGFM_PIP < .5,]
	gene_hits_df = gene_df[gene_df$TGFM_PIP >= .5,]
	gene_null_df = gene_df[gene_df$TGFM_PIP < .5,]



	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]
	red_colors=brewer.pal(n = 9, name = "Reds")
	blue_colors=brewer.pal(n = 9, name = "Blues")
	#red_colors[3]
	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.0) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=2.7) + 
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.8) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=3.0) + 
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  figure_theme() +
		  labs(x=paste0("Position [MB]\nChromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name) + 
		  theme(legend.position="bottom")

	return(pp)
}

make_tgfm_manhatten_plot_for_given_window_poster_ready <- function(snp_df, gene_df, trait_name, train_name_readable, window_name, gene_x_nudge=0.4, variant_x_nudge=-0.4) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	gene_df$tissue_name <- recode(gene_df$tissue_name, Cells_EBV_transformed_lymphocytes="Lymphocytes",Whole_Blood="Whole Blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	gene_df$labeler = paste0(gene_df$gene_id, ":", as.character(gene_df$tissue_name), "\n(PIP: ", round(gene_df$TGFM_PIP, digits=2),")")
	#maxy = max(c(max(gene_df$gwas_neg_log10_p_mean),max(snp_df$gwas_neg_log10_p)))
	gene_maxy= max(gene_df$gwas_neg_log10_p_mean)
	variant_maxy= max(snp_df$gwas_neg_log10_p)

	snp_df$labeler = paste0(snp_df$rs_id, "\n(PIP: ", round(snp_df$TGFM_PIP, digits=2),")")


	# Split both snps and genes into hits and not
	snp_hits_df = snp_df[snp_df$TGFM_PIP >= .5,]
	snp_null_df = snp_df[snp_df$TGFM_PIP < .5,]
	gene_hits_df = gene_df[gene_df$TGFM_PIP >= .5,]
	gene_null_df = gene_df[gene_df$TGFM_PIP < .5,]



	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]
	red_colors=brewer.pal(n = 9, name = "Reds")
	blue_colors=brewer.pal(n = 9, name = "Blues")
	#red_colors[3]
	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.5) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=4.7) +
		  #geom_text_repel(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = variant_x_nudge, nudge_y=variant_maxy*.1) +
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=2.3) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=5.0) + 
		  geom_text_repel(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=3.8, nudge_x = gene_x_nudge, nudge_y=gene_maxy*.1) +
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  poster_theme() +
		  labs(x=paste0("Position [MB] Chromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name_readable) + 
		  theme(legend.position="bottom")

	return(pp)
}



make_tgfm_manhatten_plot_for_given_window_paper_ready <- function(snp_df, gene_df, trait_name, train_name_readable, window_name, gene_x_nudge=0.4, variant_x_nudge=-0.4) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	gene_df$tissue_name <- recode(gene_df$tissue_name, Cells_EBV_transformed_lymphocytes="Lymphocytes",Whole_Blood="Whole Blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="Skin_Sun",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="Lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	gene_df$labeler = paste0(gene_df$gene_id, ":", as.character(gene_df$tissue_name), "\n(PIP: ", round(gene_df$TGFM_PIP, digits=2),")")
	#maxy = max(c(max(gene_df$gwas_neg_log10_p_mean),max(snp_df$gwas_neg_log10_p)))
	gene_maxy= max(gene_df$gwas_neg_log10_p_mean)
	variant_maxy= max(snp_df$gwas_neg_log10_p)

	snp_df$labeler = paste0(snp_df$rs_id, "\n(PIP: ", round(snp_df$TGFM_PIP, digits=2),")")


	# Split both snps and genes into hits and not
	snp_hits_df = snp_df[snp_df$TGFM_PIP >= .5,]
	snp_null_df = snp_df[snp_df$TGFM_PIP < .5,]
	gene_hits_df = gene_df[gene_df$TGFM_PIP >= .5,]
	gene_null_df = gene_df[gene_df$TGFM_PIP < .5,]



	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]
	red_colors=brewer.pal(n = 9, name = "Reds")
	blue_colors=brewer.pal(n = 9, name = "Blues")
	#red_colors[3]
	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.0) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=2.7) +
		  #geom_text_repel(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = variant_x_nudge, nudge_y=variant_maxy*.1) +
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.8) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=3.0) + 
		  geom_text_repel(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = gene_x_nudge, nudge_y=gene_maxy*.1) +
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  figure_theme() +
		  labs(x=paste0("Position [MB] Chromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name_readable) + 
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



if (FALSE) {
# Loop through specific examples
n_examples = dim(specific_examples_data)[1]
for (example_iter in 1:n_examples) {
	# Extract relevent info for this example
	trait_name <- specific_examples_data$trait[example_iter]
	trait_name_readable <- specific_examples_data$trait[example_iter]
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
	tgfm_manhatten_plot <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)
	output_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_tgfm_manhattan.pdf")
	ggsave(tgfm_manhatten_plot, file=output_file, width=7.2, height=4.0, units="in")
}
}



if (FALSE) {
######################
# Example 1
# ENSG00000163599_T8	ENSG00000163599	CTLA4	T8	0.8723686818022242	0.8909160835269218	2:202010553:205010553
trait_name = "disease_AID_ALL"
trait_name_readable = "Autoimmune"
window_name="2:202010553:205010553"
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
tgfm_manhatten_plot1 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)


######################
# Example 2
#ENSG00000130592_ncM	ENSG00000130592	LSP1	ncM	0.6625057303168707	0.9948638260669503	11:70863:3070863
trait_name = "blood_MONOCYTE_COUNT"
trait_name_readable = "Monocyte count"
window_name="11:70863:3070863"
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
tgfm_manhatten_plot2 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)



######################
# Example 3
#ENSG00000163421_cM	ENSG00000163421	PROK2	cM	0.5626783205874496	0.57881536884434	3:70018518:73018518
trait_name = "blood_MEAN_PLATELET_VOL"
trait_name_readable = "Mean platelet volume (MPV)"
window_name="3:70018518:73018518"
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
tgfm_manhatten_plot3 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)


######################
# Example 3
#ENSG00000163421_cM	ENSG00000163421	PROK2	cM	0.5626783205874496	0.57881536884434	3:70018518:73018518
trait_name = "blood_MEAN_PLATELET_VOL"
trait_name_readable = "Mean platelet volume (MPV)"
window_name="3:70018518:73018518"
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
tgfm_manhatten_plot3 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, gene_x_nudge=-.4)

######################
# Example 4
#ENSG00000189403_B	ENSG00000189403	HMGB1	B	0.6274887573736335	0.6274887573736335	13:29445954:32445954
trait_name = "blood_MEAN_CORPUSCULAR_HEMOGLOBIN"
trait_name_readable = "Mean corpuscular hemoglobin (MCH)"
window_name="13:29445954:32445954"
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
tgfm_manhatten_plot4 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)


####################
# Joint plot
legender = get_legend(tgfm_manhatten_plot1)
joint_manhatten <- plot_grid(plot_grid(tgfm_manhatten_plot1 + theme(legend.position="none"), tgfm_manhatten_plot2+ theme(legend.position="none"),tgfm_manhatten_plot3+ theme(legend.position="none"), tgfm_manhatten_plot4+ theme(legend.position="none"), ncol=2, labels=c("a","b","c","d")), legender, ncol=1, rel_heights=c(1, .14))
output_file <- paste0(visualize_specific_tgfm_examples_dir, "figure7.pdf")
ggsave(joint_manhatten, file=output_file, width=7.2, height=5.0, units="in")
}














if (FALSE) {
######################
# Example 1
# blood_RED_COUNT ENSG00000235169.7_Whole_Blood   SMIM1   ENSG00000235169.7       Whole_Blood     1:2010176:5010176
trait_name <- "blood_RED_COUNT"
trait_name_readable <- "Red blood cell count"
window_name <- "1:2010176:5010176"
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
tgfm_manhatten_plot1 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, variant_x_nudge=-.5)



######################
# Example 2
#disease_HYPOTHYROIDISM_SELF_REP ENSG00000115705.20_Thyroid      TPO     ENSG00000115705.20      Thyroid 2:10553:3010553
trait_name <- "disease_HYPOTHYROIDISM_SELF_REP"
trait_name_readable <- "Hypothyroidism"
window_name <- "2:10553:3010553"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot2 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)


######################
# Example 3
#lung_FEV1FVCzSMOKE      ENSG00000118257.16_Lung NRP2    ENSG00000118257.16      Lung    2:204010553:207010553
trait_name <- "lung_FEV1FVCzSMOKE"
trait_name_readable <- "FEV1:FVC"
window_name <- "2:204010553:207010553"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")
#link_df <- read.table(link_df_input_file, header=TRUE, sep="\t")
link_df=0

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot3 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)


######################
# Example 4
#disease_T2D     ENSG00000180176.14_Pancreas     TH      ENSG00000180176.14      Pancreas        11:1070863:4070863
trait_name="disease_T2D"
trait_name_readable="Type 2 Diabetes"
window_name="11:1070863:4070863"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")


# Make manhatten plot for a given window/trait
tgfm_manhatten_plot4 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, gene_x_nudge=-.4, variant_x_nudge=.8)


####################
# Joint plot
legender = get_legend(tgfm_manhatten_plot1)
joint_manhatten <- plot_grid(plot_grid(tgfm_manhatten_plot1 + theme(legend.position="none"), tgfm_manhatten_plot2+ theme(legend.position="none"), tgfm_manhatten_plot3+ theme(legend.position="none"),tgfm_manhatten_plot4+ theme(legend.position="none"), ncol=2, labels=c("a","b","c","d")), legender, ncol=1, rel_heights=c(1, .09))
output_file <- paste0(visualize_specific_tgfm_examples_dir, "four_example_tgfm_manhattan.pdf")
ggsave(joint_manhatten, file=output_file, width=7.2, height=5.2, units="in")



######################
# Example 5
#disease_HYPERTENSION_DIAGNOSED  ENSG00000159640.15_Adrenal_Gland        ACE     ENSG00000159640.15      Adrenal_Gland   17:62150508:65150508
trait_name <- "disease_HYPERTENSION_DIAGNOSED"
trait_name_readable <- "Hypertension"
window_name <- "17:62150508:65150508"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")


# Make manhatten plot for a given window/trait
tgfm_manhatten_plot5 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)


########################
# Example 6
# disease_ALLERGY_ECZEMA_DIAGNOSED        ENSG00000172818.9_Cells_EBV-transformed_lymphocytes     OVOL1   ENSG00000172818.9       Cells_EBV-transformed_lymphocytes       11:64070863:67070863
trait_name <- "disease_ALLERGY_ECZEMA_DIAGNOSED"
trait_name_readable <- "Eczema"
window_name <- "11:64070863:67070863"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot6 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, variant_x_nudge=.4)

########################
# Example 7
# disease_HYPERTENSION_DIAGNOSED  ENSG00000182263.13_Artery_Heart FIGN    ENSG00000182263.13      Artery_Heart    2:162010553:165010553
trait_name <- "disease_HYPERTENSION_DIAGNOSED"
trait_name_readable <- "Hypertension"
window_name <- "2:162010553:165010553"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot7 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, gene_x_nudge=-.4,variant_x_nudge=.4)






####################
# Joint plot
legender = get_legend(tgfm_manhatten_plot1)
joint_manhatten <- plot_grid(plot_grid(tgfm_manhatten_plot1 + theme(legend.position="none"), tgfm_manhatten_plot2+ theme(legend.position="none"), tgfm_manhatten_plot3+ theme(legend.position="none"),tgfm_manhatten_plot4+ theme(legend.position="none"), tgfm_manhatten_plot6+ theme(legend.position="none"), tgfm_manhatten_plot5+ theme(legend.position="none"), ncol=2, labels=c("a","b","c","d", "e", "f")), legender, ncol=1, rel_heights=c(1, .09))
output_file <- paste0(visualize_specific_tgfm_examples_dir, "six_example_tgfm_manhattan.pdf")
ggsave(joint_manhatten, file=output_file, width=7.2, height=6.5, units="in")


####################
# Joint plot
legender = get_legend(tgfm_manhatten_plot1)
joint_manhatten <- plot_grid(plot_grid(tgfm_manhatten_plot1 + theme(legend.position="none"), tgfm_manhatten_plot2+ theme(legend.position="none"),tgfm_manhatten_plot4+ theme(legend.position="none"), tgfm_manhatten_plot6+ theme(legend.position="none"), tgfm_manhatten_plot5+ theme(legend.position="none"), tgfm_manhatten_plot7+ theme(legend.position="none"), ncol=2, labels=c("a","b","c","d", "e", "f")), legender, ncol=1, rel_heights=c(1, .09))
output_file <- paste0(visualize_specific_tgfm_examples_dir, "six_example_tgfm_manhattan_v2.pdf")
ggsave(joint_manhatten, file=output_file, width=7.2, height=6.5, units="in")

}



########################
# Example 6
# disease_ALLERGY_ECZEMA_DIAGNOSED        ENSG00000172818.9_Cells_EBV-transformed_lymphocytes     OVOL1   ENSG00000172818.9       Cells_EBV-transformed_lymphocytes       11:64070863:67070863
trait_name <- "disease_ALLERGY_ECZEMA_DIAGNOSED"
trait_name_readable <- "Eczema"
window_name <- "11:64070863:67070863"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_posterplot <- make_tgfm_manhatten_plot_for_given_window_poster_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name)
output_file <- paste0(visualize_specific_tgfm_examples_dir, "eczema_poster_example.pdf")
ggsave(tgfm_manhatten_posterplot, file=output_file, width=7.2, height=6.0, units="in")




