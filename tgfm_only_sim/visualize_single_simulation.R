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

make_power_med_h2_se_barplot <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=power)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Power") +
  		theme(legend.position="bottom") 
  	return(p)
}

make_type_1_error_med_h2_se_barplot <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=type_1_error)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=type_1_error_lb, ymax=type_1_error_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="FDR") +
  		theme(legend.position="bottom")  +
  		geom_hline(yintercept=.05, linetype=2)

  	return(p)
}

 make_avg_h2_se_barplot <- function(df) {

 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=total_h2)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=total_h2_lb, ymax=total_h2_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Total h2") +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=.3, linetype=2)
  	return(p)
}


make_avg_fraction_h2_se_barplot <- function(df) {
	 df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=frac_h2)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=frac_h2_lb, ymax=frac_h2_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Fraction h2 mediated\nby gene expression") +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=.1, linetype=2)
  	return(p)
}

make_avg_per_tissue_h2_se_barplot_for_single_class <- function(df, titler, sim_value, min_height, max_height) {
	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	df$tissue_number = factor(df$tissue_number)
	p<-ggplot(data=df, aes(x=tissue_number, y=med_h2, fill=eqtl_sample_size)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Tissue number", y="h2 mediated by gene expression", fill="eQTL sample size", title=titler) +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=sim_value, linetype=2) +
  		ylim(min_height, max_height) +
  		theme(plot.title = element_text(hjust = 0.5))
  	return(p)



}

make_avg_per_tissue_h2_se_barplot <- function(df) {

	causal_df <- df[as.character(df$causal_status)=="causal",]
	null_df <- df[as.character(df$causal_status)=="null",]

	min_height = min(df$med_h2_lb)
	max_height = max(df$med_h2_ub)
	max_height = max(c(max_height, .0152))



	causal_plot <- make_avg_per_tissue_h2_se_barplot_for_single_class(causal_df, "Causal tissues", .015, min_height, max_height)
	null_plot <- make_avg_per_tissue_h2_se_barplot_for_single_class(null_df, "Null tissues", .0, min_height, max_height)
	legender <- get_legend(causal_plot)

	joint_plot <- plot_grid(causal_plot +theme(legend.position="none"), null_plot +theme(legend.position="none")+labs(y=""), ncol=2, rel_widths=c(.29,.8))

	joint_plot2 <- plot_grid(joint_plot, legender, ncol=1, rel_heights=c(1,.2))
	return(joint_plot2)
}


make_n_detected_genes_se_barplot <- function(df) {
	eqtl_sample_sizes <- c(100, 200, 300, 500, 1000)

	eqtl_ss_arr <- c()
	gene_type_arr <- c()
	fraction_arr <- c()
	fraction_lb_arr <- c()
	fraction_ub_arr <- c()

	for (eqtl_ss_iter in 1:length(eqtl_sample_sizes)) {
		eqtl_ss <- eqtl_sample_sizes[eqtl_ss_iter]
		small_df <- df[df$eQTL_sample_size==eqtl_ss,]
		n_detected_h <- sum(small_df$n_detected_heritable_genes)
		n_h <- sum(small_df$n_heritable_genes)
		n_detected_nh <- sum(small_df$n_detected_non_heritable_genes)
		n_nh <- sum(small_df$n_non_heritable_genes)

		frac1 <- n_detected_h/n_h
		frac2 <- n_detected_nh/n_nh
		frac1_se <- sqrt((frac1)*(1.0-frac1)/(n_h))
		frac2_se <- sqrt((frac2)*(1.0-frac2)/(n_nh))
		frac1_lb = frac1 - (1.96*frac1_se)
		frac1_ub = frac1 + (1.96*frac1_se)
		frac2_lb = frac2 - (1.96*frac2_se)
		frac2_ub = frac2 + (1.96*frac2_se)

		eqtl_ss_arr <- c(eqtl_ss_arr, eqtl_ss, eqtl_ss)
		gene_type_arr <- c(gene_type_arr, "heritable", "non-heritable")
		fraction_arr <- c(fraction_arr, frac1, frac2)
		fraction_lb_arr <- c(fraction_lb_arr, frac1_lb, frac2_lb)
		fraction_ub_arr <- c(fraction_ub_arr, frac1_ub, frac2_ub)

	}
	df2 <- data.frame(eQTL_sample_size=factor(eqtl_ss_arr, levels=eqtl_sample_sizes), gene_type=gene_type_arr, fraction=fraction_arr, frac_lb=fraction_lb_arr, frac_ub=fraction_ub_arr)

	p<-ggplot(data=df2, aes(x=eQTL_sample_size, y=fraction, fill=gene_type)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=frac_lb, ymax=frac_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="Fraction of genes detected", fill="") 
  	return(p)
}

make_tgfm_pip_fdr_plot_varying_eqtl_sample_and_element_class <- function(df2, pip_threshold, ln_pi_method) {
	df <- df2[(df2$genetic_element_class)!="all",]
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb

	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(100,200,300,500,1000, "Inf"))

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0(ln_pi_method, " prior / PIP=", pip_threshold))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12))
  	return(p)
}

make_tgfm_pip_fdr_plot_varying_prior_method_and_element_class <- function(df, pip_threshold, ln_pi_methods, eQTL_sample_size) {
	df <- df[(df$genetic_element_class)!="all",]
	df <- df[(df$genetic_element_class)!="all",]
	df = df[as.character(df$ln_pi_method) %in% ln_pi_methods,]
	df = df[df$eQTL_sample_size == eQTL_sample_size,]

	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb
	df$ln_pi_method = factor(df$ln_pi_method, levels=ln_pi_methods)

	p<-ggplot(data=df, aes(x=ln_pi_method, y=fdr, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Prior method", y="FDR", fill="", title=paste0("eqtl sample size: ", eQTL_sample_size))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)
}

make_tgfm_pip_fdr_plot_varying_detected_gene_and_eqtl_sample_size <- function(calibration_df, calibration_df2, pip_threshold, ln_pi_method, initialization_version) {
	eqtl_ss_arr <- c()
	n_element_arr <- c()
	fdr_arr <- c()
	fdr_lb_arr <- c()
	fdr_ub_arr <- c()
	eval_version_arr <- c()

	eqtl_ss_arr <- c(calibration_df$eQTL_sample_size, calibration_df2$eQTL_sample_size)
	n_element_arr <- c(calibration_df$n_elements, calibration_df2$n_elements)
	fdr_arr <- c(1.0 - calibration_df$coverage, 1.0-calibration_df2$coverage)
	fdr_lb_arr <- c(1.0 - calibration_df$coverage_ub, 1.0 - calibration_df2$coverage_ub)
	fdr_ub_arr <- c(1.0 - calibration_df$coverage_lb, 1.0 - calibration_df2$coverage_lb)
	eval_version_arr <- c(rep("standard", length(calibration_df$eQTL_sample_size)), rep("causal_genes_detected", length(calibration_df2$eQTL_sample_size)))

	df <- data.frame(eQTL_sample_size=factor(eqtl_ss_arr), n_elements=n_element_arr, fdr=fdr_arr, fdr_lb=fdr_lb_arr, fdr_ub=fdr_ub_arr, eval_version=factor(eval_version_arr, levels=c("standard", "causal_genes_detected")))
	
	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=eval_version)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP: ", pip_threshold, " / Prior: ", ln_pi_method))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)
  }


 make_tgfm_pip_fdr_plot_varying_twas_method_and_element_class <- function(df, pip_threshold, ln_pi_method, eqtl_sample_size, initialization_version) {
 	df <- df[(df$genetic_element_class)!="all",]
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb
	df$twas_method = factor(df$twas_method, levels=c("susie_pmces", "susie_distr", "fusion_lasso_pmces"))

	p<-ggplot(data=df, aes(x=twas_method, y=fdr, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="", y="FDR", fill="", title=paste0("eqtl sample size: ", eqtl_sample_size))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)			
 }


make_tgfm_pip_fdr_plot_varying_initialization_and_element_class <- function(df, pip_threshold, ln_pi_method, eqtl_sample_size) {
	df <- df[(df$genetic_element_class)!="all",]
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb
	df$initialization_version = factor(df$initialization_version, levels=c("null", "variant_only", "best"))

	p<-ggplot(data=df, aes(x=initialization_version, y=fdr, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="TGFM Initialization", y="FDR", fill="", title=paste0("eqtl sample size: ", eqtl_sample_size, " / PIP: ", pip_threshold, " / Prior: ", ln_pi_method))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)		
}

make_tgfm_pip_fdr_plot_varying_sample_size_and_n_causal_elements <- function(df, pip_threshold) {
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb

	df$eQTL_sample_size = factor(df$eQTL_sample_size)
	df$n_causal_bin = factor(df$n_causal_bin, levels=c("1_8", "9_15", "16_30"))


	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=n_causal_bin)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="# causal elements", title=paste0("PIP: ", pip_threshold))  +
  		geom_hline(yintercept=(1.0-as.numeric(pip_threshold)), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)		
}

make_tgfm_pip_power_plot_varying_twas_method_and_element_class <- function(df, pip_threshold, ln_pi_method, eqtl_sample_size) {
	df$twas_method = factor(df$twas_method, levels=c("susie_pmces", "susie_distr", "fusion_lasso_pmces"))
	p<-ggplot(data=df, aes(x=twas_method, y=power, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="", y="Power", fill="", title=paste0("eqtl sample size: ", eqtl_sample_size))  +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


  	return(p)
}

make_tgfm_pip_power_plot_varying_initialization_and_element_class <- function(df, pip_threshold, ln_pi_method, eqtl_sample_size) {
	df$initialization_version = factor(df$initialization_version, levels=c("null", "variant_only", "best"))
	p<-ggplot(data=df, aes(x=initialization_version, y=power, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="TGFM Initialization", y="Power", fill="", title=paste0("eqtl sample size: ", eqtl_sample_size, " / PIP: ", pip_threshold, " / Prior: ", ln_pi_method))  +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


  	return(p)

}


make_tgfm_pip_power_plot_varying_prior_method_and_element_class <- function(df, pip_threshold, ln_pi_methods, eQTL_sample_size) {
	df = df[as.character(df$ln_pi_method) %in% ln_pi_methods,]
	df = df[df$eQTL_sample_size == eQTL_sample_size,]

	df$ln_pi_method = factor(df$ln_pi_method, levels=ln_pi_methods)


	p<-ggplot(data=df, aes(x=ln_pi_method, y=power, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Prior method", y="Power", fill="", title=paste0("eqtl sample size: ", eQTL_sample_size))  +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


  	return(p)

}


make_tgfm_pip_power_plot_varying_eqtl_sample_and_element_class <- function(df, pip_threshold, ln_pi_method) {
	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(100,200,300,500,1000, "Inf"))

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0(ln_pi_method, " prior / PIP=", pip_threshold))  +
  		theme(plot.title = element_text(hjust = 0.5,size=12))

  	return(p)
}

make_number_of_elements_per_window_histogram <- function(df, ele_class) {
	df = df[as.character(df$genetic_element_class) == ele_class,]

	pp <- ggplot(df, aes(x=n_causal_genetic_elements))+
  		geom_histogram(color="darkblue", fill="lightblue", binwidth = 1.0) +
  		labs(x=paste0("# Causal genetic ", ele_class, "s / genome window"), title=ele_class) +
  		figure_theme()
  	return(pp)

}

make_number_of_elements_per_window_fp_strat_histogram <- function(df, pip_threshold) {
	info <- wilcox.test(df$n_causal_genetic_elements[df$fp_gene=="True"], df$n_causal_genetic_elements[df$fp_gene=="False"], alternative = "two.sided")
	pvalue = info$p.value
	pp <- ggplot(df, aes(x=n_causal_genetic_elements, color=fp_gene)) +
  		geom_histogram(aes(y = ..density..),fill="white", alpha=0.5, position="identity", binwidth=1.0) +
  		figure_theme() +
  		labs(x=paste0("# Causal genetic elements / genome window"), color="Window has\nFP gene",title=paste0("PIP: ", pip_threshold, " / p-value=", pvalue)) 
  	return(pp)
}

make_gene_causal_element_ld_score_fp_strat_histogram <- function(df) {
	info <- wilcox.test(df$causal_element_ld_score[df$false_positive_gene=="True"], df$causal_element_ld_score[df$false_positive_gene=="False"], alternative = "two.sided")
	pvalue = info$p.value
	pp <- ggplot(df, aes(x=causal_element_ld_score, color=false_positive_gene)) +
  		geom_histogram(aes(y = ..density..),fill="white", alpha=0.5, position="identity") +
  		figure_theme() +
  		labs(x=paste0("Causal genetic element LD score / gene"), color="FP gene",title=paste0("PIP: 0.9", " / p-value=", pvalue)) 
  	return(pp)
}


make_tgfm_pip_fdr_plot_varying_n_causal_elements_and_twas_method <- function(df, pip_threshold, eqtl_sample_size) {
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb


	df$n_causal_bin = factor(df$n_causal_bin, levels=c("1_8", "9_15", "16_30"))
	df$twas_method = factor(df$twas_method, levels=c("susie_pmces", "bootstrapped_susie_distr"))

	p<-ggplot(data=df, aes(x=n_causal_bin, y=fdr, fill=twas_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="# causal elements", y="FDR", fill="", title=paste0("PIP: ", pip_threshold, " / eqtl ss: ", eqtl_sample_size))  +
  		geom_hline(yintercept=(1.0-as.numeric(pip_threshold)), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)		
}


#######################
# Command line args
#######################
global_simulation_name_string = args[1]
simulated_organized_results_dir = args[2]
visualize_simulated_results_dir = args[3]





#####################################################################
# TGFM
#####################################################################

#####################################################################
# Calibration at varying pip thresholds while just showing gene FDR stratefied by sample size and n_causal-genetic elements
#####################################################################
ln_pi_method="uniform"
initialization_version="best"
resid_var_method="False"
pip_threshold_arr <- c("0.5", "0.9", "0.95")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration_by_n_causal.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
	calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
	calibration_df = calibration_df[as.character(calibration_df$residual_variance) == resid_var_method,]
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) == "gene",]


	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_sample_size_and_n_causal_elements(calibration_df, pip_threshold)
	plots[[pip_iter]] = calibration_barplot
}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_gene_fdr_vary_sample_size_and_n_causal_elements.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")




#####################################################################
# Calibration at varying pip thresholds while just showing gene FDR stratefied by method type and n_causal-genetic elements
#####################################################################
ln_pi_method="uniform"
initialization_version="best"
resid_var_method="False"
pip_threshold_arr <- c("0.3", "0.5", "0.9", "0.95")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_fixed_500_method_compare__tgfm_pip_", pip_threshold, "_calibration_by_n_causal.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
	calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
	calibration_df = calibration_df[as.character(calibration_df$residual_variance) == resid_var_method,]
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) == "gene",]
	calibration_df = calibration_df[as.character(calibration_df$twas_method) != "marginal_distr",]
	calibration_df = calibration_df[as.character(calibration_df$twas_method) != "bootstrapped_marginal_distr",]

	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_n_causal_elements_and_twas_method(calibration_df, pip_threshold, "500") + theme(legend.position="bottom")
	plots[[pip_iter]] = calibration_barplot
}
joint <- plot_grid(plotlist=plots, ncol=2)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_gene_fdr_vary_twas_method_and_n_causal_elements.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")






if (FALSE) {
#####################################################################
# General sim data
#####################################################################
# Plot causal-genetic element LD score/ gene stratefied by false positive status
gene_causal_element_ld_score_file = paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_tgfm_inf_susie_pmces_gene_causal_genetic_element_ld_scores.txt")
gene_causal_element_ldscore_df <- read.table(gene_causal_element_ld_score_file, header=TRUE)
gene_causal_element_ld_score_fp_strat_histogram <- make_gene_causal_element_ld_score_fp_strat_histogram(gene_causal_element_ldscore_df)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_causal_element_ld_score_fp_strat_histogram.pdf")

ggsave(gene_causal_element_ld_score_fp_strat_histogram, file=output_file, width=7.2, height=3.7, units="in")





# Plot distribution of number of causal elements per window
organized_number_causal_elements_per_region_file = paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_number_of_causal_genetic_elements_per_region.txt")
n_causal_ele_df <- read.table(organized_number_causal_elements_per_region_file, header=TRUE)

ele_class <- "variant"
n_var_per_window_histogram <- make_number_of_elements_per_window_histogram(n_causal_ele_df, ele_class)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_n_causal_", ele_class, "_per_window_histogram.pdf")
ggsave(n_var_per_window_histogram, file=output_file, width=7.2, height=2.7, units="in")

ele_class <- "gene"
n_var_per_window_histogram <- make_number_of_elements_per_window_histogram(n_causal_ele_df, ele_class)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_n_causal_", ele_class, "_per_window_histogram.pdf")
ggsave(n_var_per_window_histogram, file=output_file, width=7.2, height=2.7, units="in")


# Plot distribution of number of causal elements stratefied by false positive windows
pip_threshold = "0.9"
n_causal_elements_by_false_positive_file = paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_tgfm_inf_eqtl_pip_",pip_threshold,"_causal_genetic_density_vs_gene_false_positives.txt")
n_causal_ele_fp_df <- read.table(n_causal_elements_by_false_positive_file, header=TRUE)
n_var_per_window_fp_strat_histogram <- make_number_of_elements_per_window_fp_strat_histogram(n_causal_ele_fp_df, pip_threshold)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_n_causal_elements_per_window_fp_strat_", pip_threshold, "_histogram.pdf")
ggsave(n_var_per_window_fp_strat_histogram, file=output_file, width=7.2, height=2.7, units="in")
pip_threshold = "0.95"
n_causal_elements_by_false_positive_file = paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_tgfm_inf_eqtl_pip_",pip_threshold,"_causal_genetic_density_vs_gene_false_positives.txt")
n_causal_ele_fp_df <- read.table(n_causal_elements_by_false_positive_file, header=TRUE)
n_var_per_window_fp_strat_histogram <- make_number_of_elements_per_window_fp_strat_histogram(n_causal_ele_fp_df, pip_threshold)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_n_causal_elements_per_window_fp_strat_", pip_threshold, "_histogram.pdf")
ggsave(n_var_per_window_fp_strat_histogram, file=output_file, width=7.2, height=2.7, units="in")
pip_threshold = "0.99"
n_causal_elements_by_false_positive_file = paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_tgfm_inf_eqtl_pip_",pip_threshold,"_causal_genetic_density_vs_gene_false_positives.txt")
n_causal_ele_fp_df <- read.table(n_causal_elements_by_false_positive_file, header=TRUE)
n_var_per_window_fp_strat_histogram <- make_number_of_elements_per_window_fp_strat_histogram(n_causal_ele_fp_df, pip_threshold)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_n_causal_elements_per_window_fp_strat_", pip_threshold, "_histogram.pdf")
ggsave(n_var_per_window_fp_strat_histogram, file=output_file, width=7.2, height=2.7, units="in")

}


if (FALSE) {
#####################################################################
# TGFM
#####################################################################

#####################################################################
# Calibration at PIP threshold of .9 at varying twas versions for variant and gene assuming single prior at fixed sample size
#####################################################################
ln_pi_method="uniform"
eqtl_ss_arr <- c("100", "300", "500", "1000")
plots <- list()
for (eqtl_ss_iter in 1:length(eqtl_ss_arr)) {
initialization_version="best"
pip_threshold=0.9
resid_var_method="False"
eqtl_sample_size = eqtl_ss_arr[eqtl_ss_iter]
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df = calibration_df[as.character(calibration_df$eQTL_sample_size) == eqtl_sample_size,]
calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
calibration_df = calibration_df[as.character(calibration_df$residual_variance) == resid_var_method,]

# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_twas_method_and_element_class(calibration_df, pip_threshold, ln_pi_method, eqtl_sample_size, initialization_version)
plots[[eqtl_ss_iter]] = calibration_barplot
}
joint <- plot_grid(plotlist=plots, ncol=2)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_fdr_vary_twas_method_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")


#####################################################################
# Power at PIP threshold of .9 at twas versions for variant and gene assuming single prior and fixed sample size
#####################################################################

pip_threshold=0.9
ln_pi_method="uniform"
initialization_version="best"
resid_var_method="False"

plots <- list()

eqtl_ss_arr <- c("100", "300", "500", "1000")
plots <- list()
for (eqtl_ss_iter in 1:length(eqtl_ss_arr)) {
# Load in data
eqtl_sample_size = eqtl_ss_arr[eqtl_ss_iter]

power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
power_df <- read.table(power_file, header=TRUE)
power_df = power_df[as.character(power_df$ln_pi_method) == ln_pi_method,]
power_df = power_df[as.character(power_df$eQTL_sample_size) == eqtl_sample_size,]
power_df = power_df[as.character(power_df$initialization_version) == initialization_version,]
power_df = power_df[as.character(power_df$residual_variance) == resid_var_method,]



# Make plot
power_barplot <- make_tgfm_pip_power_plot_varying_twas_method_and_element_class(power_df, pip_threshold, ln_pi_method, eqtl_sample_size)
plots[[eqtl_ss_iter]] = power_barplot
}
joint <- plot_grid(plotlist=plots, ncol=2)

# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_power_vary_twas_method_with_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")


}





if (FALSE) {


#####################################################################
# Calibration at PIP threshold of .9 at varying twas versions for variant and gene assuming single prior at fixed sample size
# Calibration for genes not gene-tissue pairs
#####################################################################
ln_pi_method="uniform"
eqtl_ss_arr <- c("100", "300", "500", "1000")
plots <- list()
for (eqtl_ss_iter in 1:length(eqtl_ss_arr)) {
initialization_version="best"
pip_threshold=0.9
eqtl_sample_size = eqtl_ss_arr[eqtl_ss_iter]
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_only_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df = calibration_df[as.character(calibration_df$eQTL_sample_size) == eqtl_sample_size,]
calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_twas_method_and_element_class(calibration_df, pip_threshold, ln_pi_method, eqtl_sample_size, initialization_version)
plots[[eqtl_ss_iter]] = calibration_barplot
}
joint <- plot_grid(plotlist=plots, ncol=2)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_gene_only_fdr_vary_twas_method_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")








#####################################################################
# Calibration at PIP threshold of .9 at varying eQTL sample sizes for variant and gene assuming single prior
#####################################################################
pip_thresholds = c(0.5, 0.9)
for (pip_threshold_iter in 1:length(pip_thresholds)) {
pip_threshold=pip_thresholds[pip_threshold_iter]
ln_pi_method="uniform"
initialization_version="best"
twas_method="fusion_lasso_pmces"
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
calibration_df = calibration_df[as.character(calibration_df$twas_method) == twas_method,]

# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_and_element_class(calibration_df, pip_threshold, ln_pi_method)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_fdr_with_", ln_pi_method, "_prior_", pip_threshold, "_pip_", twas_method,".pdf")
ggsave(calibration_barplot, file=output_file, width=7.2, height=4.5, units="in")

}

#####################################################################
# Calibration at PIP threshold of .9 at varying eQTL sample sizes for variant and gene assuming single prior
#####################################################################
pip_thresholds = c(0.5, 0.9)
for (pip_threshold_iter in 1:length(pip_thresholds)) {
pip_threshold=pip_thresholds[pip_threshold_iter]
ln_pi_method="uniform"
initialization_version="best"
twas_method="fusion_lasso_pmces"
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_only_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
calibration_df = calibration_df[as.character(calibration_df$twas_method) == twas_method,]

# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_and_element_class(calibration_df, pip_threshold, ln_pi_method)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_gene_only_fdr_with_", ln_pi_method, "_prior_", pip_threshold, "_pip_", twas_method,".pdf")
ggsave(calibration_barplot, file=output_file, width=7.2, height=4.5, units="in")

}




}























if (FALSE) {
#####################################################################
# Calibration at PIP threshold of .9 at varying eQTL sample sizes for variant and gene assuming single prior
#####################################################################
pip_thresholds = c(0.5, 0.9)
for (pip_threshold_iter in 1:length(pip_thresholds)) {
pip_threshold=pip_thresholds[pip_threshold_iter]
ln_pi_method="uniform"
initialization_version="best"
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_and_element_class(calibration_df, pip_threshold, ln_pi_method)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_fdr_with_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(calibration_barplot, file=output_file, width=7.2, height=4.5, units="in")
}

#####################################################################
# Power at PIP threshold of .9 at varying eQTL sample sizes for variant and gene assuming single prior
#####################################################################
pip_thresholds = c(0.5, 0.9)
for (pip_threshold_iter in 1:length(pip_thresholds)) {
pip_threshold=pip_thresholds[pip_threshold_iter]
ln_pi_method="uniform"
initialization_version="best"
# Load in data
power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
power_df <- read.table(power_file, header=TRUE)
power_df = power_df[as.character(power_df$ln_pi_method) == ln_pi_method,]
power_df = power_df[as.character(power_df$initialization_version) == initialization_version,]
# Make plot
power_barplot <- make_tgfm_pip_power_plot_varying_eqtl_sample_and_element_class(power_df, pip_threshold, ln_pi_method)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_power_with_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(power_barplot, file=output_file, width=7.2, height=4.5, units="in")
}


#####################################################################
# Calibration at PIP threshold of .9 at varying initialization versions for variant and gene assuming single prior at fixed sample size
#####################################################################
pip_thresholds = c(0.5, 0.9, 0.99)
for (pip_threshold_iter in 1:length(pip_thresholds)) {
pip_threshold=pip_thresholds[pip_threshold_iter]
ln_pi_method="uniform"
eqtl_sample_size="Inf"
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df = calibration_df[as.character(calibration_df$eQTL_sample_size) == eqtl_sample_size,]
# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_initialization_and_element_class(calibration_df, pip_threshold, ln_pi_method, eqtl_sample_size)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_fdr_with_", ln_pi_method, "_prior_", eqtl_sample_size, "_eQTL_ss_", pip_threshold, "_pip.pdf")
ggsave(calibration_barplot, file=output_file, width=7.2, height=3.7, units="in")
}

#####################################################################
# Power at PIP threshold of .9 at initialization versions for variant and gene assuming single prior and fixed sample size
#####################################################################
pip_thresholds = c(0.5, 0.9, 0.99)
for (pip_threshold_iter in 1:length(pip_thresholds)) {
pip_threshold=pip_thresholds[pip_threshold_iter]
ln_pi_method="uniform"
eqtl_sample_size="Inf"
# Load in data
power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
power_df <- read.table(power_file, header=TRUE)
power_df = power_df[as.character(power_df$ln_pi_method) == ln_pi_method,]
power_df = power_df[as.character(power_df$eQTL_sample_size) == eqtl_sample_size,]
# Make plot
power_barplot <- make_tgfm_pip_power_plot_varying_initialization_and_element_class(power_df, pip_threshold, ln_pi_method, eqtl_sample_size)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_power_with_", ln_pi_method, "_prior_", eqtl_sample_size, "_eQTL_ss_", pip_threshold, "_pip.pdf")
ggsave(power_barplot, file=output_file, width=7.2, height=3.7, units="in")
}

#####################################################################
# Calibration at PIP threshold of .9 at for gene while varying eQTL sample size and whether genes needs to be detected
#####################################################################
pip_threshold=0.9
ln_pi_method="uniform"
eqtl_sample_sizes=c("100", "200", "300", "500", "1000")
initialization_version="best"
genetic_element_type="gene"

calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df = calibration_df[as.character(calibration_df$eQTL_sample_size) %in% eqtl_sample_sizes,]
calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) == genetic_element_type,]
calibration_df = calibration_df[as.character(calibration_df$initialization_version) == initialization_version,]
calibration_file2 <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration_where_causal_gene_is_detected.txt")
calibration_df2 <- read.table(calibration_file2, header=TRUE)
calibration_df2 = calibration_df2[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
calibration_df2 = calibration_df2[as.character(calibration_df2$eQTL_sample_size) %in% eqtl_sample_sizes,]
calibration_df2 = calibration_df2[as.character(calibration_df2$genetic_element_class) == genetic_element_type,]
calibration_df2 = calibration_df2[as.character(calibration_df2$initialization_version) == initialization_version,]

# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_detected_gene_and_eqtl_sample_size(calibration_df, calibration_df2, pip_threshold, ln_pi_method, initialization_version)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_fdr_varying_gene_detection_evaluation_with_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(calibration_barplot, file=output_file, width=7.2, height=3.7, units="in")





}





















if (FALSE) {


#####################################################################
# Calibration at PIP threshold of .9 at fixed eqtl sample sample size varying prior and for variant and gene assuming single prior
#####################################################################
pip_threshold=0.9
eQTL_sample_size=1000
ln_pi_methods <- c("uniform", "variant_v_gene_only_1e-08", "sparse_estimate_1e-08", "point_estimate_1e-08")
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
calibration_df <- read.table(calibration_file, header=TRUE)
# Make plot
calibration_barplot <- make_tgfm_pip_fdr_plot_varying_prior_method_and_element_class(calibration_df, pip_threshold, ln_pi_methods, eQTL_sample_size)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_fdr_with_eqtl_ss_", eQTL_sample_size,".pdf")
ggsave(calibration_barplot, file=output_file, width=7.2, height=4.5, units="in")




#####################################################################
# Power at PIP threshold of .9 at fixed eqtl sample sample size varying prior and for variant and gene assuming single prior
#####################################################################
pip_threshold=0.9
eQTL_sample_size=1000
ln_pi_methods <- c("uniform", "variant_v_gene_only_1e-08", "sparse_estimate_1e-08", "point_estimate_1e-08")
# Load in data
power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
power_df <- read.table(power_file, header=TRUE)
# Make plot
power_barplot <- make_tgfm_pip_power_plot_varying_prior_method_and_element_class(power_df, pip_threshold, ln_pi_methods, eQTL_sample_size)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_power_with_eqtl_ss_", eQTL_sample_size,".pdf")
ggsave(power_barplot, file=output_file, width=7.2, height=4.5, units="in")


}







