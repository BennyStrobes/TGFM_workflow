args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(reshape2)
library(ggbeeswarm)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


mean_se_barplot_of_prop_drug_target_gene_score_binned_by_tgfm_pip <- function(df) {

	pops_mean_vec <- c()
	pops_mean_se_vec <- c()
	bin_names_vec <- c()

	threshold_lb <- 0.0
	threshold_ub <- .1
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	threshold_lb <- .1
	threshold_ub <- .25
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))


	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	threshold_lb <- .5
	threshold_ub <- .75
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	threshold_lb <- .75
	threshold_ub <- 1.0
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	df2 <- data.frame(pops=pops_mean_vec, pops_se=pops_mean_se_vec, tgfm_bin=factor(bin_names_vec))

	p <- ggplot(df2) +
    		geom_bar( aes(x=tgfm_bin, y=pops), stat="identity", fill="skyblue", alpha=0.7) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="Proportion of genes\nin drug target gene set", x="") +
    		geom_errorbar( aes(x=tgfm_bin, ymin=pops-(1.96*pops_se), ymax=pops+(1.96*pops_se)), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    		figure_theme()
    return(p)




}

mean_se_barplot_of_prop_drug_target_gene_score_binned_by_tgfm_pip_and_twas_z <- function(df) {

	pops_mean_vec <- c()
	pops_mean_se_vec <- c()
	bin_names_vec <- c()
	method_type_vec <- c()
	n_indices = 0

	thresher = 0.0

	sorted_twas = sort(df$max_abs_twas_z)



	#TGFM
	threshold_lb <- 0.0
	threshold_ub <- .1
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices+1)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	method_type_vec <- c(method_type_vec, "TWAS")

	#TGFM
	threshold_lb <- .1
	threshold_ub <- .25
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices+1)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	method_type_vec <- c(method_type_vec, "TWAS")

	#TGFM
	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices+1)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	method_type_vec <- c(method_type_vec, "TWAS")

	#TGFM
	threshold_lb <- .5
	threshold_ub <- 1.0
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	prop = mean(df$drug_target_gene[indices])
	prop_se = sqrt((prop*(1.0-prop))/length(df$drug_target_gene[indices]))
	pops_mean_vec <- c(pops_mean_vec, prop)
	pops_mean_se_vec <- c(pops_mean_se_vec, prop_se)
	method_type_vec <- c(method_type_vec, "TWAS")


	df2 <- data.frame(pops=pops_mean_vec, pops_se=pops_mean_se_vec, tgfm_bin=factor(bin_names_vec), methody=factor(method_type_vec, levels=c("TWAS", "TGFM")))

	p <- ggplot(df2, aes(x=tgfm_bin, y=pops, fill=methody)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		scale_fill_manual(values=c("grey", "plum3")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="Proportion of genes\nin drug target gene set", x="", fill="") +
    		geom_errorbar( aes(ymin=pops-(1.96*pops_se), ymax=pops+(1.96*pops_se)), width=0.2, colour="grey50", position=position_dodge(.9)) +
    		figure_theme()
    return(p)
}





#######################
# Command line args
########################
working_dir <- args[1]

# Load in summary data
summary_df <- read.table(paste0(working_dir, "drug_target_gene_set_overlap_summary_file.txt"), header=TRUE)





######################
# average and standard error of proportion of drug-target genes binned by TGFM PIP
output_file <- paste0(working_dir, "mean_se_barplot_prop_drug_target_gene_binned_by_tgfm_pip_x_trait.pdf")
barplot <- mean_se_barplot_of_prop_drug_target_gene_score_binned_by_tgfm_pip(summary_df)
ggsave(barplot, file=output_file, width=7.2, height=4.7, units="in")

######################
# average and standard error of proportion of drug-target genes binned by TGFM PIP and twas z
output_file <- paste0(working_dir, "mean_se_barplot_prop_drug_target_gene_binned_by_tgfm_pip_and_twas_z_x_trait.pdf")
barplot <- mean_se_barplot_of_prop_drug_target_gene_score_binned_by_tgfm_pip_and_twas_z(summary_df)
ggsave(barplot, file=output_file, width=7.2, height=4.7, units="in")



