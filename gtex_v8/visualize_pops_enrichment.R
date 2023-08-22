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





make_scatterplot_of_twas_z_vs_tgfm_pip_colored_by_pops_score <- function(df) {

	df$twas_log_10_pvalue = -log10(2*pnorm(q=df$max_abs_twas_z, lower.tail=FALSE))

	pops_threshold = -sort(-df$pops_score)[round(.05*length(df$pops_score))]

	df$pops_score[df$pops_score > pops_threshold] = pops_threshold
	df$pops_boolean = df$pops > pops_threshold
	p <- ggplot(df, aes(x=twas_log_10_pvalue, y=max_tgfm_pip, color=pops_score)) +
  		geom_point(size=.1) +
  		figure_theme() 
  		#labs(title=trait_name, x="TGFM expected #\ncausal genes (PIP >=.2)", y="TGLR expected mediated h2")
  	return(p)

}

make_boxplot_of_pops_score_binned_by_tgfm_pip <- function(df) {
	pops_vec <- c()
	bin_names <- c()


	threshold_lb <- 0.0
	threshold_ub <- .1
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	pops_vec <- c(pops_vec, df$pops_score[indices])
	bin_names <- c(bin_names, rep(paste0(threshold_lb," <= TGFM PIP < ", threshold_ub), length(df$pops_score[indices])))

	threshold_lb <- .1
	threshold_ub <- .25
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	pops_vec <- c(pops_vec, df$pops_score[indices])
	bin_names <- c(bin_names, rep(paste0(threshold_lb," <= TGFM PIP < ", threshold_ub), length(df$pops_score[indices])))

	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	pops_vec <- c(pops_vec, df$pops_score[indices])
	bin_names <- c(bin_names, rep(paste0(threshold_lb," <= TGFM PIP < ", threshold_ub), length(df$pops_score[indices])))

	threshold_lb <- .5
	threshold_ub <- 1.0
	indices = (df$max_tgfm_pip >= threshold_lb) & (df$max_tgfm_pip < threshold_ub)
	pops_vec <- c(pops_vec, df$pops_score[indices])
	bin_names <- c(bin_names, rep(paste0(threshold_lb," <= TGFM PIP < ", threshold_ub), length(df$pops_score[indices])))




	df2 = data.frame(pops=pops_vec, tgfm_bin=factor(bin_names))

	p <- ggplot(df2, aes(x=tgfm_bin, y=pops_vec)) + 
  		geom_boxplot() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		labs(x="", y="POPS Score") +
  		figure_theme()

  	return(p)

}

mean_se_barplot_of_pops_score_binned_by_tgfm_pip_and_twas_z <- function(df) {
	pops_mean_vec <- c()
	pops_mean_se_vec <- c()
	bin_names_vec <- c()
	method_type_vec <- c()
	n_indices = 0

	thresher = 0.0

	sorted_twas = sort(df$max_abs_twas_z)
	sorted_tgfm_twas = sort(df$tgfm_abs_gene_pmces)


	###############################
	threshold_lb <- 0.0
	threshold_ub <- .1
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	#TGFM
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices+1)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TWAS")	
	#TGFMTWAS
	#start_index = n_indices
	#n_indices = n_indices + sum(indices)
	threshold_lb = sorted_tgfm_twas[(start_index+1)]
	threshold_ub = sorted_tgfm_twas[(n_indices+1)]
	indices = (df$tgfm_abs_gene_pmces >= threshold_lb) & (df$tgfm_abs_gene_pmces < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM-Z")	

	###############################
	threshold_lb <- .1
	threshold_ub <- .25
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	method_type_vec <- c(method_type_vec, "TGFM")
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices+1)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TWAS")
	#TGFMTWAS
	#start_index = n_indices
	#n_indices = n_indices + sum(indices)
	threshold_lb = sorted_tgfm_twas[(start_index+1)]
	threshold_ub = sorted_tgfm_twas[(n_indices+1)]
	indices = (df$tgfm_abs_gene_pmces >= threshold_lb) & (df$tgfm_abs_gene_pmces < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM-Z")	

	###############################
	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	method_type_vec <- c(method_type_vec, "TGFM")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices+1)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TWAS")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)
	#TGFMTWAS
	#start_index = n_indices
	#n_indices = n_indices + sum(indices)
	threshold_lb = sorted_tgfm_twas[(start_index+1)]
	threshold_ub = sorted_tgfm_twas[(n_indices+1)]
	indices = (df$tgfm_abs_gene_pmces >= threshold_lb) & (df$tgfm_abs_gene_pmces < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM-Z")	


	###############################
	threshold_lb <- .5
	threshold_ub <- .7
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices+1)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TWAS")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)
	#TGFMTWAS
	#start_index = n_indices
	#n_indices = n_indices + sum(indices)
	threshold_lb = sorted_tgfm_twas[(start_index+1)]
	threshold_ub = sorted_tgfm_twas[(n_indices+1)]
	indices = (df$tgfm_abs_gene_pmces >= threshold_lb) & (df$tgfm_abs_gene_pmces < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM-Z")	

	###############################
	threshold_lb <- .7
	threshold_ub <- 0.9
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TWAS")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)	
	#TGFMTWAS
	#start_index = n_indices
	#n_indices = n_indices + sum(indices)
	threshold_lb = sorted_tgfm_twas[(start_index+1)]
	threshold_ub = sorted_tgfm_twas[(n_indices+1)]
	indices = (df$tgfm_abs_gene_pmces >= threshold_lb) & (df$tgfm_abs_gene_pmces < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM-Z")	


	###############################
	threshold_lb <- 0.9
	threshold_ub <- 1.0
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)
	#TWAS
	start_index = n_indices
	n_indices = n_indices + sum(indices)
	threshold_lb = sorted_twas[(start_index+1)]
	threshold_ub = sorted_twas[(n_indices)]
	indices = (df$max_abs_twas_z >= threshold_lb) & (df$max_abs_twas_z < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TWAS")
	fraction = sum(df$pops_score[indices] < thresher)/sum(indices)
	print(fraction)	
	#TGFMTWAS
	#start_index = n_indices
	#n_indices = n_indices + sum(indices)
	threshold_lb = sorted_tgfm_twas[(start_index+1)]
	threshold_ub = sorted_tgfm_twas[(n_indices)]
	indices = (df$tgfm_abs_gene_pmces >= threshold_lb) & (df$tgfm_abs_gene_pmces < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	#bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))
	method_type_vec <- c(method_type_vec, "TGFM-Z")	




	df2 <- data.frame(pops=pops_mean_vec, pops_se=pops_mean_se_vec, tgfm_bin=factor(bin_names_vec), methody=factor(method_type_vec, levels=c("TWAS","TGFM-Z", "TGFM")))


	print(df2)

	p <- ggplot(df2, aes(x=tgfm_bin, y=pops, fill=methody)) +
    		geom_bar(stat="identity", position=position_dodge()) +
    		scale_fill_manual(values=c("grey", "red","plum3")) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="Average POPS score", x="", fill="") +
    		geom_errorbar( aes(ymin=pops-(1.96*pops_se), ymax=pops+(1.96*pops_se)), width=0.2, colour="grey50", position=position_dodge(.9)) +
    		figure_theme()
    return(p)


}


mean_se_barplot_of_pops_score_binned_by_tgfm_pip <- function(df) {

	pops_mean_vec <- c()
	pops_mean_se_vec <- c()
	bin_names_vec <- c()

	threshold_lb <- 0.0
	threshold_ub <- .1
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	threshold_lb <- .1
	threshold_ub <- .25
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))


	threshold_lb <- .25
	threshold_ub <- .5
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	threshold_lb <- .5
	threshold_ub <- .7
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))

	threshold_lb <- .7
	threshold_ub <- .9
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))


	threshold_lb <- .1
	threshold_ub <- 1.0
	indices = (df$tgfm_gene_pip >= threshold_lb) & (df$tgfm_gene_pip < threshold_ub)
	pops_mean_vec <- c(pops_mean_vec, mean(df$pops_score[indices]))
	pops_mean_se_vec <- c(pops_mean_se_vec, sd(df$pops_score[indices])/sqrt(length(df$pops_score[indices])))
	bin_names_vec <- c(bin_names_vec, paste0(threshold_lb," <= TGFM PIP < ", threshold_ub))


	df2 <- data.frame(pops=pops_mean_vec, pops_se=pops_mean_se_vec, tgfm_bin=factor(bin_names_vec))

	green_color=brewer.pal(n = 9, name = "Greens")[6]

	print("HI")

	p <- ggplot(df2) +
    		geom_bar( aes(x=tgfm_bin, y=pops), stat="identity", fill=green_color, alpha=0.9) +
    		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
    		labs(y="Average POPS score", x="") +
    		geom_errorbar( aes(x=tgfm_bin, ymin=pops-(1.96*pops_se), ymax=pops+(1.96*pops_se)), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    		figure_theme()
    return(p)


}

#######################
# Command line args
########################
working_dir <- args[1]

# Load in summary data
summary_df <- read.table(paste0(working_dir, "cross_traits_pops_tgfm_enrichment_summary.txt"), header=TRUE)





######################
# Scatter plot of TGFM-PIP vs TWAS Z-score colored by POPS score
output_file <- paste0(working_dir, "scatterplot_twas_z_vs_tgfm_pip_colored_by_pops_score_x_trait.pdf")
#scatterplot <- make_scatterplot_of_twas_z_vs_tgfm_pip_colored_by_pops_score(summary_df)
#ggsave(scatterplot, file=output_file, width=7.2, height=6.7, units="in")



######################
# Boxplot of pops-score binned by TGFM PIP
output_file <- paste0(working_dir, "boxplot_pops_score_binned_by_tgfm_pip_x_trait.pdf")
#boxplot <- make_boxplot_of_pops_score_binned_by_tgfm_pip(summary_df)
#ggsave(boxplot, file=output_file, width=7.2, height=4.7, units="in")

######################
# average and standard error of mean of pops-score binned by TGFM PIP
output_file <- paste0(working_dir, "mean_se_barplot_pops_score_binned_by_tgfm_pip_x_trait.pdf")
barplot <- mean_se_barplot_of_pops_score_binned_by_tgfm_pip(summary_df)
ggsave(barplot, file=output_file, width=7.2, height=4.7, units="in")

######################
# average and standard error of mean of pops-score binned by TGFM PIP and twas-z
output_file <- paste0(working_dir, "mean_se_barplot_pops_score_binned_by_tgfm_pip_and_twas_z_x_trait.pdf")
barplot <- mean_se_barplot_of_pops_score_binned_by_tgfm_pip_and_twas_z(summary_df)
ggsave(barplot, file=output_file, width=7.2, height=4.7, units="in")


