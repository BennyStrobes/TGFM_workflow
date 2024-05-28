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

get_tgfm_manhatten_plot_for_given_window_paper_data <- function(snp_df, gene_df, trait_name, train_name_readable, window_name) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	#gene_df$tissue_name <- recode(gene_df$tissue_name, Thyroid="thyroid",Artery_Aorta="artery aorta", Brain_Cerebellum="brain cerebellum",Cells_EBV_transformed_lymphocytes="lymphocytes",Whole_Blood="whole blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="skin (sun exposed)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	#gene_df$tissue_name = recode(gene_df$tissue_name, T4="CD4", T8="CD8")

	position_vec <- c(gene_df$gene_tss, snp_df$snp_position)
	name_vec <- c(as.character(gene_df$gene_name), as.character(snp_df$snp_name))
	association_vec <- c(gene_df$gwas_neg_log10_p_mean, snp_df$gwas_neg_log10_p)
	pip_vec <- c(gene_df$TGFM_PIP, snp_df$TGFM_PIP)
	trait_vec <- c(rep(trait_name, length(gene_df$gene_tss)), rep(trait_name, length(snp_df$snp_position)))
	window_vec <- c(rep(window_name, length(gene_df$gene_tss)), rep(window_name, length(snp_df$snp_position)))
	class_vec <- c(rep("Gene-tissue", length(gene_df$gene_tss)), rep("Variant", length(snp_df$snp_position)))


	df <- data.frame(trait=trait_vec, genomic_window=window_vec, genetic_element_class=class_vec, genetic_element=name_vec, association_neg_log10_p=association_vec, TGFM_PIP=pip_vec, position=position_vec)

	return(df)
}

make_tgfm_manhatten_plot_for_given_window_paper_ready <- function(snp_df, gene_df, trait_name, train_name_readable, window_name, twas_thresh=4.2e-7, gene_x_nudge=0.4, variant_x_nudge=-0.4, add_nominal_sig_lines=FALSE, add_twas_nominal_sig_lines=FALSE, add_nm_variant_text=FALSE, variant_scaler=.1, gene_scaler=.1, gene_y_nudge=0.0) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	gene_df$tissue_name <- recode(gene_df$tissue_name, Thyroid="thyroid",Artery_Aorta="artery aorta", Brain_Cerebellum="brain cerebellum",Cells_EBV_transformed_lymphocytes="lymphocytes",Whole_Blood="whole blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="skin (sun exposed)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord", T4="CD4", T8="CD8")

	# DIAGNOSTIC STUFF
	print(trait_name)
	best_gene_tissue_index = which.max(gene_df$TGFM_PIP)
	best_gene_tissue_name = gene_df$gene_id_tissue_name[best_gene_tissue_index]
	best_gene_name = gene_df
	best_gene_tissue_pip = gene_df$TGFM_PIP[best_gene_tissue_index]
	best_gene_tissue_twas_p = 10^(-gene_df$gwas_neg_log10_p_mean[best_gene_tissue_index])
	print(paste0(best_gene_tissue_name, " : ", best_gene_tissue_pip, " : ", best_gene_tissue_twas_p))

	gene_df2 = gene_df[abs(gene_df$gene_tss_mb - gene_df$gene_tss_mb[best_gene_tissue_index]) < 1,]
	twas_p = 10^(-gene_df2$gwas_neg_log10_p_mean)
	print(dim(gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print((gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print(snp_df[snp_df$TGFM_PIP >= .5,])

	print((gene_df2[(twas_p <= twas_thresh),]))


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

	print(twas_thresh)

	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]
	red_colors=brewer.pal(n = 9, name = "Reds")
	blue_colors=brewer.pal(n = 9, name = "Blues")
	#red_colors[3]
	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.0) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=3.0) 

	if (add_nm_variant_text) {
		pp = pp + geom_text_repel(data=snp_hits_df, size=2.0, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP), nudge_x = variant_x_nudge, nudge_y=variant_maxy*variant_scaler)
	}
	pp = pp + scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.8) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=3.0) + 
		  geom_text_repel(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = gene_x_nudge, nudge_y=gene_maxy*gene_scaler + gene_y_nudge) +
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  figure_theme() +
		  labs(x=paste0("Position [MB] Chromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name_readable) + 
		  theme(legend.position="bottom")
	if (add_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')+
		geom_hline(yintercept=-log10(5*10^-8), color=brewer.pal(n = 9, name = "Blues")[4], linetype='dotted')
	}
	if (add_twas_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')
	}


	return(pp)
}


make_tgfm_manhatten_plot_for_given_window_paper_ready_custom_SBP_AO2 <- function(snp_df, gene_df, trait_name, train_name_readable, window_name, twas_thresh=4.2e-7, gene_x_nudge=0.4, variant_x_nudge=-0.4, add_nominal_sig_lines=FALSE, add_twas_nominal_sig_lines=FALSE, add_nm_variant_text=FALSE, variant_scaler=.1, gene_scaler=.1, gene_y_nudge=0.0) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	gene_df$tissue_name <- recode(gene_df$tissue_name, Thyroid="thyroid",Artery_Aorta="artery aorta", Brain_Cerebellum="brain cerebellum",Cells_EBV_transformed_lymphocytes="lymphocytes",Whole_Blood="whole blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="skin (sun exposed)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord", T4="CD4", T8="CD8")

	# DIAGNOSTIC STUFF
	print(trait_name)
	best_gene_tissue_index = which.max(gene_df$TGFM_PIP)
	best_gene_tissue_name = gene_df$gene_id_tissue_name[best_gene_tissue_index]
	best_gene_name = gene_df
	best_gene_tissue_pip = gene_df$TGFM_PIP[best_gene_tissue_index]
	best_gene_tissue_twas_p = 10^(-gene_df$gwas_neg_log10_p_mean[best_gene_tissue_index])
	print(paste0(best_gene_tissue_name, " : ", best_gene_tissue_pip, " : ", best_gene_tissue_twas_p))

	gene_df2 = gene_df[abs(gene_df$gene_tss_mb - gene_df$gene_tss_mb[best_gene_tissue_index]) < 1,]
	twas_p = 10^(-gene_df2$gwas_neg_log10_p_mean)
	print(dim(gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print((gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print(snp_df[snp_df$TGFM_PIP >= .5,])

	print((gene_df2[(twas_p <= twas_thresh),]))


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

	print(twas_thresh)

	snp_hits_df1 <- snp_hits_df[as.character(snp_hits_df$rs_id) == "rs145778816",]
	snp_hits_df2 <- snp_hits_df[as.character(snp_hits_df$rs_id) == "rs138682554",]

	gene_hits_df1 <- gene_hits_df[as.character(gene_hits_df$gene_id) == "IDH2", ]
	gene_hits_df2 <- gene_hits_df[as.character(gene_hits_df$gene_id) == "FES", ]


	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]
	red_colors=brewer.pal(n = 9, name = "Reds")
	blue_colors=brewer.pal(n = 9, name = "Blues")
	#red_colors[3]
	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.0) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=3.0) 

	if (add_nm_variant_text) {
		pp = pp + geom_text_repel(data=snp_hits_df1, size=3.0, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP), nudge_x = 0.01, nudge_y=variant_maxy*variant_scaler)
		pp = pp + geom_text_repel(data=snp_hits_df2, size=3.0, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP), nudge_x = 0.1, nudge_y=variant_maxy*variant_scaler*.07)

	}
	pp = pp + scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.8) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=3.0) + 
		  geom_text_repel(data=gene_hits_df1, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=3.0, nudge_x = 0.0, nudge_y=gene_maxy*gene_scaler + 5.8,segment.color = NA) +
		  geom_text_repel(data=gene_hits_df2, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=3.0, nudge_x = 0.16, nudge_y=gene_maxy*gene_scaler + 2.0) +
  

		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  figure_theme() +
		  labs(x=paste0("Position [MB] Chromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name_readable) + 
		  theme(legend.position="bottom")
	if (add_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')+
		geom_hline(yintercept=-log10(5*10^-8), color=brewer.pal(n = 9, name = "Blues")[4], linetype='dotted')
	}
	if (add_twas_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')
	}


	return(pp)
}



make_tgfm_manhatten_plot_for_given_window_paper_ready_custom_SBP_AO <- function(snp_df, gene_df, trait_name, train_name_readable, window_name, twas_thresh=4.2e-7, gene_x_nudge=0.4, variant_x_nudge=-0.4, add_nominal_sig_lines=FALSE, add_twas_nominal_sig_lines=FALSE, add_nm_variant_text=FALSE, variant_scaler=.1, gene_scaler=.1, gene_y_nudge=0.0) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	gene_df$tissue_name <- recode(gene_df$tissue_name, Thyroid="thyroid",Artery_Aorta="artery aorta", Brain_Cerebellum="brain cerebellum",Cells_EBV_transformed_lymphocytes="lymphocytes",Whole_Blood="whole blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="skin (sun exposed)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord", T4="CD4", T8="CD8")

	# DIAGNOSTIC STUFF
	print(trait_name)
	best_gene_tissue_index = which.max(gene_df$TGFM_PIP)
	best_gene_tissue_name = gene_df$gene_id_tissue_name[best_gene_tissue_index]
	best_gene_name = gene_df
	best_gene_tissue_pip = gene_df$TGFM_PIP[best_gene_tissue_index]
	best_gene_tissue_twas_p = 10^(-gene_df$gwas_neg_log10_p_mean[best_gene_tissue_index])
	print(paste0(best_gene_tissue_name, " : ", best_gene_tissue_pip, " : ", best_gene_tissue_twas_p))

	gene_df2 = gene_df[abs(gene_df$gene_tss_mb - gene_df$gene_tss_mb[best_gene_tissue_index]) < 1,]
	twas_p = 10^(-gene_df2$gwas_neg_log10_p_mean)
	print(dim(gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print((gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print(snp_df[snp_df$TGFM_PIP >= .5,])

	print((gene_df2[(twas_p <= twas_thresh),]))


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

	print(twas_thresh)

	snp_hits_df1 <- snp_hits_df[as.character(snp_hits_df$rs_id) == "rs145778816",]
	snp_hits_df2 <- snp_hits_df[as.character(snp_hits_df$rs_id) == "rs138682554",]

	gene_hits_df1 <- gene_hits_df[as.character(gene_hits_df$gene_id) == "IDH2", ]
	gene_hits_df2 <- gene_hits_df[as.character(gene_hits_df$gene_id) == "FES", ]


	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]
	red_colors=brewer.pal(n = 9, name = "Reds")
	blue_colors=brewer.pal(n = 9, name = "Blues")
	#red_colors[3]
	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.0) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=3.0) 

	if (add_nm_variant_text) {
		pp = pp + geom_text_repel(data=snp_hits_df1, size=2.0, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP), nudge_x = 0.01, nudge_y=variant_maxy*variant_scaler)
		pp = pp + geom_text_repel(data=snp_hits_df2, size=2.0, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP), nudge_x = 0.1, nudge_y=variant_maxy*variant_scaler*.07)

	}
	pp = pp + scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.8) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=3.0) + 
		  geom_text_repel(data=gene_hits_df1, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = -.35, nudge_y=gene_maxy*gene_scaler + 5.8,segment.color = NA) +
		  geom_text_repel(data=gene_hits_df2, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = .55, nudge_y=gene_maxy*gene_scaler + 2.0) +
  

		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  figure_theme() +
		  labs(x=paste0("Position [MB] Chromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name_readable) + 
		  theme(legend.position="bottom")
	if (add_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')+
		geom_hline(yintercept=-log10(5*10^-8), color=brewer.pal(n = 9, name = "Blues")[4], linetype='dotted')
	}
	if (add_twas_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')
	}


	return(pp)
}





make_tgfm_manhatten_plot_for_given_window_paper_ready_tmp <- function(snp_df, gene_df, trait_name, train_name_readable, window_name, twas_thresh=4.2e-7, gene_x_nudge=0.4, variant_x_nudge=-0.4, add_nominal_sig_lines=FALSE, add_twas_nominal_sig_lines=FALSE, add_nm_variant_text=FALSE, variant_scaler=.1, gene_scaler=.1, gene_y_nudge=0.0, gene_x_nudge2=.4, gene_y_nudge2=0.0) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	gene_df$tissue_name <- recode(gene_df$tissue_name, Thyroid="thyroid",Artery_Aorta="artery aorta", Brain_Cerebellum="brain cerebellum",Cells_EBV_transformed_lymphocytes="lymphocytes",Whole_Blood="whole blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="skin (sun exposed)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord", T4="CD4", T8="CD8")

	# DIAGNOSTIC STUFF
	print(trait_name)
	best_gene_tissue_index = which.max(gene_df$TGFM_PIP)
	best_gene_tissue_name = gene_df$gene_id_tissue_name[best_gene_tissue_index]
	best_gene_name = gene_df
	best_gene_tissue_pip = gene_df$TGFM_PIP[best_gene_tissue_index]
	best_gene_tissue_twas_p = 10^(-gene_df$gwas_neg_log10_p_mean[best_gene_tissue_index])
	print(paste0(best_gene_tissue_name, " : ", best_gene_tissue_pip, " : ", best_gene_tissue_twas_p))

	gene_df2 = gene_df[abs(gene_df$gene_tss_mb - gene_df$gene_tss_mb[best_gene_tissue_index]) < 1,]
	twas_p = 10^(-gene_df2$gwas_neg_log10_p_mean)
	print(dim(gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print((gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print(snp_df[snp_df$TGFM_PIP >= .5,])

	print((gene_df2[(twas_p <= twas_thresh),]))


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

	gene_hits_df1 = gene_hits_df[as.character(gene_hits_df$gene_id) == 'IDH2',]
	gene_hits_df2 = gene_hits_df[as.character(gene_hits_df$gene_id) != 'IDH2',]

	# Get name of chromosome
	chromosome_num = strsplit(as.character(window_name),split=":")[[1]][1]
	red_colors=brewer.pal(n = 9, name = "Reds")
	blue_colors=brewer.pal(n = 9, name = "Blues")
	#red_colors[3]
	pp <- ggplot() + 
		  geom_point(data=snp_null_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=1.0) + 
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=3.0) 

	if (add_nm_variant_text) {
		pp = pp + geom_text_repel(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = variant_x_nudge, nudge_y=variant_maxy*variant_scaler)
	}
	pp = pp + scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.8) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=3.0) + 
		  geom_text_repel(data=gene_hits_df1, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = gene_x_nudge, nudge_y=gene_maxy*gene_scaler + gene_y_nudge) +
		  geom_text_repel(data=gene_hits_df2, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = gene_x_nudge2, nudge_y=gene_maxy*gene_scaler + gene_y_nudge2) +

		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  figure_theme() +
		  labs(x=paste0("Position [MB] Chromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name_readable) + 
		  theme(legend.position="bottom")
	if (add_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')+
		geom_hline(yintercept=-log10(5*10^-8), color=brewer.pal(n = 9, name = "Blues")[4], linetype='dotted')
	}
	if (add_twas_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')
	}


	return(pp)
}

make_tgfm_manhatten_plot_for_given_window_ppt_ready <- function(snp_df, gene_df, trait_name, train_name_readable, window_name, twas_thresh=4.2e-7, gene_x_nudge=0.4, variant_x_nudge=-0.4, add_nominal_sig_lines=FALSE, add_twas_nominal_sig_lines=FALSE, add_nm_variant_text=TRUE, variant_scaler=.1, gene_scaler=.1) {
	# Convert snp_df and gene_df to MB
	snp_df$snp_position_mb = snp_df$snp_position/1000000.0
	gene_df$gene_tss_mb = gene_df$gene_tss/1000000.0

	gene_df$tissue_name = str_replace_all(as.character(gene_df$tissue_name), "-", "_")

	gene_df$tissue_name <- recode(gene_df$tissue_name, Thyroid="thyroid",Artery_Aorta="artery aorta", Brain_Cerebellum="brain cerebellum",Cells_EBV_transformed_lymphocytes="lymphocytes",Whole_Blood="whole blood", Adrenal_Gland="Adrenal Gland",Adipose_Subcutaneous="Adipose_Sub", Adipose_Visceral_Omentum="Adipose_Visceral", Breast_Mammary_Tissue="Breast_Mammary", Cells_Cultured_fibroblasts="Fibroblast",Heart_Atrial_Appendage="Heart_Atrial",Skin_Sun_Exposed_Lower_leg="skin (sun exposed)",Skin_Not_Sun_Exposed_Suprapubic="Skin_No_Sun", Small_Intestine_Terminal_Ileum="Small_Intestine", Brain_Anterior_cingulate_cortex_BA24="Brain_anterior_cortex", Brain_Nucleus_accumbens_basal_ganglia="Brain_basal_ganglia", Esophagus_Gastroesophageal_Junction="Esophagus_gastro_jxn", Cells_EBV_transformed_lymphocytes="lymphocytes", Brain_Spinal_cord_cervical_c_1="Brain_Spinal_cord")

	# DIAGNOSTIC STUFF
	if (FALSE) {
	print(trait_name)
	best_gene_tissue_index = which.max(gene_df$TGFM_PIP)
	best_gene_tissue_name = gene_df$gene_id_tissue_name[best_gene_tissue_index]
	best_gene_name = gene_df
	best_gene_tissue_pip = gene_df$TGFM_PIP[best_gene_tissue_index]
	best_gene_tissue_twas_p = 10^(-gene_df$gwas_neg_log10_p_mean[best_gene_tissue_index])
	print(paste0(best_gene_tissue_name, " : ", best_gene_tissue_pip, " : ", best_gene_tissue_twas_p))

	gene_df2 = gene_df[abs(gene_df$gene_tss_mb - gene_df$gene_tss_mb[best_gene_tissue_index]) < 1,]
	twas_p = 10^(-gene_df2$gwas_neg_log10_p_mean)
	print(dim(gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print((gene_df2[(twas_p <= twas_thresh) & gene_df2$TGFM_PIP <= .01,]))
	print(snp_df[snp_df$TGFM_PIP >= .5,])

	print((gene_df2[(twas_p <= twas_thresh),]))
	}


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
		  geom_point(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, color=TGFM_PIP), shape=16, size=3.0) 

	if (add_nm_variant_text) {
		#pp = pp + geom_text_repel(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = variant_x_nudge, nudge_y=variant_maxy*variant_scaler)
		pp = pp + geom_text_repel(data=snp_hits_df, aes(x=snp_position_mb, y=gwas_neg_log10_p, label=labeler, color=TGFM_PIP),nudge_x=.09,size=3.0, nudge_y=.07)
	}
	pp = pp + scale_colour_gradient(low = blue_colors[3], high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  labs(color="Variant PIP") +
		  new_scale_color() + 
		  geom_point(data=gene_null_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17,size=1.8) + 
		  geom_point(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, color=TGFM_PIP), shape=17, size=3.0) + 
		  #geom_text_repel(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=2.8, nudge_x = gene_x_nudge, nudge_y=gene_maxy*gene_scaler) +
		  geom_text_repel(data=gene_hits_df, aes(x=gene_tss_mb, y=gwas_neg_log10_p_mean, label=labeler, color=TGFM_PIP),size=3.0, nudge_x=.04, nudge_y=.4) +
		  #scale_colour_gradient(low = "lightskyblue1", high = "dodgerblue4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  scale_colour_gradient(low = red_colors[3], high = "red4", limits=c(0.0,1.0), breaks=c(0.0,.25,.5,.75,1.0),labels=c("0.0",".25",".5",".75","1.0")) +
		  figure_theme() +
		  labs(x=paste0("Position [MB] Chromosome ", chromosome_num), y=bquote(-log[10](pvalue)), color="Gene-Tissue PIP", title=trait_name_readable) + 
		  theme(legend.position="bottom")
	if (add_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')+
		geom_hline(yintercept=-log10(5*10^-8), color=brewer.pal(n = 9, name = "Blues")[4], linetype='dotted')
	}
	if (add_twas_nominal_sig_lines==TRUE) {
		pp = pp +geom_hline(yintercept=-log10(twas_thresh), color=brewer.pal(n = 9, name = "Reds")[4], linetype='dotted')
	}


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

if (FALSE) {
twas_thresh=4.2e-7
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
	tgfm_manhatten_plot <- make_tgfm_manhatten_plot_for_given_window_ppt_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE)
	output_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_tgfm_manhattan.pdf")
	ggsave(tgfm_manhatten_plot, file=output_file, width=7.2, height=3.6, units="in")
}
}
if (FALSE) {
twas_thresh=4.1e-7
######################
# Example 1
trait_name = "blood_MONOCYTE_COUNT"
trait_name_readable = "Monocyte count"
window_name="1:25010176:28010176"
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
tgfm_manhatten_plot1 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name,twas_thresh=twas_thresh, add_nominal_sig_lines=TRUE)
tgfm_manhatten_data1 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)


######################
# Example 2
trait_name = "blood_LYMPHOCYTE_COUNT"
trait_name_readable = "Lymphocyte count"
window_name="15:29795210:32795210"
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
tgfm_manhatten_plot2 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, gene_scaler=.01)
tgfm_manhatten_data2 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)



######################
# Example 3
#ENSG00000189403_B	ENSG00000189403	HMGB1	B	0.6274887573736335	0.6274887573736335	13:29445954:32445954
trait_name = "blood_MEAN_CORPUSCULAR_HEMOGLOBIN"
trait_name_readable = "Corp. hemoglobin"
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
tgfm_manhatten_plot3 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE)
tgfm_manhatten_data3 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)




######################
# Example 4
# ENSG00000163599_T8	ENSG00000163599	CTLA4	T8	0.8723686818022242	0.8909160835269218	2:202010553:205010553
trait_name = "disease_AID_ALL"
trait_name_readable = "All autoimmune"
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
tgfm_manhatten_plot4 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name,twas_thresh=twas_thresh, add_nominal_sig_lines=TRUE)
tgfm_manhatten_data4 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)




####################
# Joint plot
legender = get_legend(tgfm_manhatten_plot1)
joint_manhatten <- plot_grid(plot_grid(tgfm_manhatten_plot1 + theme(legend.position="none"), tgfm_manhatten_plot2+ theme(legend.position="none"), tgfm_manhatten_plot4+ theme(legend.position="none"),tgfm_manhatten_plot3+ theme(legend.position="none"), ncol=2, labels=c("a","b","c","d")), legender, ncol=1, rel_heights=c(1, .14))
output_file <- paste0(visualize_specific_tgfm_examples_dir, "figure8.pdf")
ggsave(joint_manhatten, file=output_file, width=7.2, height=5.0, units="in")

###################
# Merge data
merged_tgfm_manhatten_data = rbind(tgfm_manhatten_data1, tgfm_manhatten_data2, tgfm_manhatten_data3, tgfm_manhatten_data4)
supp_table_file = paste0(visualize_specific_tgfm_examples_dir, "suppTable_figure8_numerical.txt")
write.table(merged_tgfm_manhatten_data, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
print(supp_table_file)
}

if (FALSE) {
# FIGURE 6 (revision)
twas_thresh=4.2e-7
########################
# Example 1
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
tgfm_manhatten_plot1 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,variant_x_nudge=.4, add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE)
tgfm_manhatten_data1 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)



######################
# Example 2
trait_name <- "biochemistry_VitaminD"
trait_name_readable <- "Vitamin D"
window_name <- "1:16010176:19010176"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot2 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE)
tgfm_manhatten_data2 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)



######################
# Example 3
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
tgfm_manhatten_plot3 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE)
tgfm_manhatten_data3 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)


######################
# Example 4
trait_name <- "repro_MENARCHE_AGE"
trait_name_readable <- "Menarche age"
window_name <- "17:43150508:46150508"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")


# Make manhatten plot for a given window/trait
tgfm_manhatten_plot4 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_scaler=.01,gene_x_nudge=-.4)
tgfm_manhatten_data4 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)



######################
# Example 5
trait_name <- "bp_SYSTOLICadjMEDz"
trait_name_readable <- "Systolic blood pressure"
window_name <- "8:41061773:44061773"
window_name <- "15:88795210:91795210"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot6 <- make_tgfm_manhatten_plot_for_given_window_paper_ready_custom_SBP_AO(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_x_nudge=0.1, gene_x_nudge=-.5, gene_y_nudge=4.5)
tgfm_manhatten_data6 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)

######################
# Example 6
trait_name <- "bp_SYSTOLICadjMEDz"
trait_name_readable <- "Systolic blood pressure"
window_name <- "8:41061773:44061773"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot5 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE)
tgfm_manhatten_data5 <- get_tgfm_manhatten_plot_for_given_window_paper_data(snp_df, gene_df, trait_name, trait_name_readable, window_name)




####################
# Joint plot
legender = get_legend(tgfm_manhatten_plot1)
joint_manhatten <- plot_grid(plot_grid(tgfm_manhatten_plot3+ theme(legend.position="none"), tgfm_manhatten_plot1 + theme(legend.position="none"), tgfm_manhatten_plot2+ theme(legend.position="none"),tgfm_manhatten_plot6+ theme(legend.position="none"), tgfm_manhatten_plot5+ theme(legend.position="none"), tgfm_manhatten_plot4+ theme(legend.position="none"), ncol=2, labels=c("a","b","c","d", "e", "f")), legender, ncol=1, rel_heights=c(1, .09))
output_file <- paste0(visualize_specific_tgfm_examples_dir, "figure6.pdf")
print(output_file)
ggsave(joint_manhatten, file=output_file, width=7.2, height=6.5, units="in")


###################
# Merge data
merged_tgfm_manhatten_data = rbind(tgfm_manhatten_data1, tgfm_manhatten_data2, tgfm_manhatten_data3, tgfm_manhatten_data4, tgfm_manhatten_data5, tgfm_manhatten_data6)
supp_table_file = paste0(visualize_specific_tgfm_examples_dir, "suppTable_figure6.txt")
write.table(merged_tgfm_manhatten_data, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
print(supp_table_file)
}




####################
# Specifically plot fig 6d
######################
if (FALSE) {
# Example 5
twas_thresh=4.2e-7

trait_name <- "bp_SYSTOLICadjMEDz"
trait_name_readable <- "Systolic blood pressure"
window_name <- "8:41061773:44061773"
window_name <- "15:88795210:91795210"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot6 <- make_tgfm_manhatten_plot_for_given_window_paper_ready_custom_SBP_AO2(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_x_nudge=0.0, gene_x_nudge=0.0, gene_y_nudge=4.5)
output_file <- paste0(visualize_specific_tgfm_examples_dir, "figure6d.pdf")
ggsave(tgfm_manhatten_plot6, file=output_file, width=7.2, height=3.6, units="in")
print(output_file)
}

#################
# Supp figure
if (FALSE) {
twas_thresh=4.2e-7

######################
# Example 1
# blood_RED_COUNT ENSG00000235169.7_Whole_Blood   SMIM1   ENSG00000235169.7       Whole_Blood     1:2010176:5010176
trait_name <- "blood_RED_COUNT"
trait_name_readable <- "Red blood cell count"
window_name <- "1:2010176:5010176"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot1 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_scaler=.01)

######################
# Example 2
# ENSG00000205978.5_Liver	ENSG00000205978	NYNRIN	Liver	0.7019168867675932	0.9772922945911715	14:23223582:26223582
trait_name <- "biochemistry_Cholesterol"
trait_name_readable <- "Cholesterol"
window_name <- "14:23223582:26223582"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot2 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_scaler=.01)


######################
# Example 3
# ENSG00000159640.15_Adrenal_Gland	ENSG00000159640	ACE	Adrenal_Gland	0.607149800587197	0.8157302610278647	17:62150508:65150508
if (FALSE) {
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
tgfm_manhatten_plot3 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_scaler=.01)
}
trait_name <- "biochemistry_VitaminD"
trait_name_readable <- "Vitamin D"
window_name <- "15:56795210:59795210"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot3 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_scaler=.01)

######################
# Example 4
# ENSG00000118257.16_Lung	ENSG00000118257	NRP2	Lung	0.97931106037432	0.9947268843600326	2:204010553:207010553
trait_name <- "lung_FEV1FVCzSMOKE"
trait_name_readable <- "FEV1:FVC"
window_name <- "2:204010553:207010553"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot4 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_scaler=.01)



######################
# Example 5
# ENSG00000166228.8_Pancreas	ENSG00000166228	PCBD1	Pancreas	0.6100154395324857	0.7599280466673501	10:69014743:72014743
# ENSG00000138083.4_Pancreas	ENSG00000138083	SIX3	Pancreas	0.6697373643679546	0.6697373643679546	2:43010553:46010553
trait_name <- "biochemistry_HbA1c"
trait_name_readable <- "HbA1c"
window_name <- "2:43010553:46010553"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
tgfm_manhatten_plot5 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, add_nm_variant_text=TRUE, variant_scaler=.01)



######################
# Example 6
trait_name <- "bp_SYSTOLICadjMEDz"
trait_name_readable <- "Systolic blood pressure"
window_name <- "18:59011148:62011148"
# Input files
snp_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_snp_df.txt")
gene_df_input_file <- paste0(visualize_specific_tgfm_examples_dir, trait_name, "_", window_name, "_gene_df.txt")
# Load in input
snp_df <- read.table(snp_df_input_file, header=TRUE, sep="\t")
gene_df <- read.table(gene_df_input_file, header=TRUE, sep="\t")

# Make manhatten plot for a given window/trait
#tgfm_manhatten_plot6 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_twas_nominal_sig_lines=TRUE)
tgfm_manhatten_plot6 <- make_tgfm_manhatten_plot_for_given_window_paper_ready(snp_df, gene_df, trait_name, trait_name_readable, window_name, twas_thresh=twas_thresh,add_nominal_sig_lines=TRUE, gene_scaler=.07)



####################
# Joint plot
legender = get_legend(tgfm_manhatten_plot1)
joint_manhatten <- plot_grid(plot_grid(tgfm_manhatten_plot3+ theme(legend.position="none"),tgfm_manhatten_plot6+ theme(legend.position="none"),tgfm_manhatten_plot1+ theme(legend.position="none"), tgfm_manhatten_plot2 + theme(legend.position="none"), tgfm_manhatten_plot4+ theme(legend.position="none"), tgfm_manhatten_plot5+ theme(legend.position="none"), ncol=2, labels=c("a","b","c","d", "e", "f")), legender, ncol=1, rel_heights=c(1, .09))
output_file <- paste0(visualize_specific_tgfm_examples_dir, "supp_figure_six_example_tgfm_manhattan.pdf")
print(output_file)
ggsave(joint_manhatten, file=output_file, width=7.2, height=6.5, units="in")
}

