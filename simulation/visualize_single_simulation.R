args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(RColorBrewer)
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

make_power_med_h2_se_barplot_across_thresholds <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
 	df$threshold = factor(df$threshold)
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=power, fill=threshold)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Power") +
  		theme(legend.position="bottom") 
  	return(p)
}


make_power_med_h2_se_barplot_gaussian_approximation <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
 	#df = df[df$threshold==threshold_val,]
 	#df$threshold = factor(df$threshold)
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=power)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Power")   	
  	return(p)
}

make_power_med_h2_se_barplot_at_single_threshold <- function(df, threshold_val) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
 	df = df[df$threshold==threshold_val,]
 	#df$threshold = factor(df$threshold)
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=power)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Power")   	
  	return(p)
}



make_type_1_error_med_h2_se_barplot <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=type_1_error)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=type_1_error_lb, ymax=type_1_error_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Type 1 Error") +
  		theme(legend.position="bottom")  +
  		geom_hline(yintercept=.05, linetype=2)

  	return(p)
}

make_type_1_error_med_h2_se_barplot_across_thresholds <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
 	df$threshold = factor(df$threshold)
	p<-ggplot(data=df, aes(x=eqtl_sample_size, fill=threshold, y=type_1_error)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=type_1_error_lb, ymax=type_1_error_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Type 1 Error") +
  		theme(legend.position="bottom")  +
  		geom_hline(yintercept=.05, linetype=2)

  	return(p)
}
make_type_1_error_med_h2_se_barplot_at_single_threshold <- function(df, threshold_val) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
 	#df$threshold = factor(df$threshold)
 	df = df[df$threshold==threshold_val,]
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=type_1_error)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=type_1_error_lb, ymax=type_1_error_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Type 1 Error") +
  		theme(legend.position="bottom")  +
  		geom_hline(yintercept=.05, linetype=2)

  	return(p)
}

make_type_1_error_med_h2_se_barplot_gaussian_approximation <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
 	#df$threshold = factor(df$threshold)
 	#df = df[df$threshold==threshold_val,]
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=type_1_error)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=type_1_error_lb, ymax=type_1_error_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Type 1 Error") +
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

make_avg_per_tissue_tau_se_barplot_for_single_class <- function(df, titler, min_height, max_height) {
	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	df$tissue_number = factor(df$tissue_number)
	p<-ggplot(data=df, aes(x=tissue_number, y=tau, fill=eqtl_sample_size)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=tau_lb, ymax=tau_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Tissue number", y="TGLR tau", fill="eQTL sample size", title=titler) +
  		theme(legend.position="bottom") +
  		#geom_hline(yintercept=sim_value, linetype=2) +
  		ylim(min_height, max_height) +
  		theme(plot.title = element_text(hjust = 0.5))
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
make_avg_per_tissue_fraction_causal_se_barplot_for_single_class <- function(df, titler, min_height, max_height) {
	#df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(300,500,1000))
	df$tissue_number = factor(df$tissue_number)
	p<-ggplot(data=df, aes(x=tissue_number, y=fraction_causal, fill=eqtl_sample_size)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fraction_causal_lb, ymax=fraction_causal_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Tissue number", y="Fraction causal", fill="eQTL sample size", title=titler) +
  		theme(legend.position="bottom") +
  		#geom_hline(yintercept=sim_value, linetype=2) +
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

make_avg_per_tissue_tau_se_barplot <- function(df) {

	causal_df <- df[as.character(df$causal_status)=="causal",]
	null_df <- df[as.character(df$causal_status)=="null",]

	min_height = min(df$tau_lb)
	min_height = min(c(min_height, 0.0))
	max_height = max(df$tau_ub)
	#max_height = max(c(max_height, .0152))



	causal_plot <- make_avg_per_tissue_tau_se_barplot_for_single_class(causal_df, "Causal tissues", min_height, max_height)
	null_plot <- make_avg_per_tissue_tau_se_barplot_for_single_class(null_df, "Null tissues", min_height, max_height)
	legender <- get_legend(causal_plot)

	joint_plot <- plot_grid(causal_plot +theme(legend.position="none"), null_plot +theme(legend.position="none")+labs(y=""), ncol=2, rel_widths=c(.29,.8))

	joint_plot2 <- plot_grid(joint_plot, legender, ncol=1, rel_heights=c(1,.2))
	return(joint_plot2)
}

make_avg_per_tissue_fraction_causal_se_barplot <- function(df) {
	causal_df <- df[as.character(df$causal_status)=="causal",]
	null_df <- df[as.character(df$causal_status)=="null",]

	min_height = min(df$fraction_causal_lb)
	min_height = min(c(min_height, 0.0))
	max_height = max(df$fraction_causal_ub)
	#max_height = max(c(max_height, .0152))



	causal_plot <- make_avg_per_tissue_fraction_causal_se_barplot_for_single_class(causal_df, "Causal tissues", min_height, max_height)
	null_plot <- make_avg_per_tissue_fraction_causal_se_barplot_for_single_class(null_df, "Null tissues", min_height, max_height)
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
	df$twas_method = factor(df$twas_method, levels=c("susie_pmces", "susie_distr"))

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

make_tgfm_pip_power_plot_varying_twas_method_and_element_class <- function(df, pip_threshold, ln_pi_method, eqtl_sample_size) {
	df$twas_method = factor(df$twas_method, levels=c("susie_pmces", "susie_distr"))
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


make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_twas_method <- function(df, pip_threshold, ln_pi_method, element_class) {
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb

	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(300,500,1000))
	
	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=twas_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		 scale_fill_manual(values=c("#F8766D", "#56B4E9"))+

  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP=", pip_threshold, " / ", ln_pi_method, " / ", element_class))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12))
  	return(p)
}

make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_element_class <- function(df, pip_threshold, ln_pi_method, twas_method) {
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb

	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(300,500,1000))
	
	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=genetic_element_class)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("#F8766D", "#999999"))+

  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP=", pip_threshold, " / ", ln_pi_method, " / ", twas_method))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12))
  	return(p)
}

make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_prior_method <- function(df, pip_threshold, new_method, element_class) {
	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb

	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(300,500,1000))

	df$ln_pi_method = factor(df$ln_pi_method, levels=c("uniform", "tglr_bootstrapped_nonnegative_pmces", "tglr_bootstrapped_nonnegative_sampler", "pmces_uniform_iterative_variant_gene_prior_pip_level_pmces", "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"))

	df=df[!is.na(df$ln_pi_method),]	
	df=df[!is.na(df$eQTL_sample_size),]	
	df$ln_pi_method = recode(df$ln_pi_method, uniform="Uniform", tglr_bootstrapped_nonnegative_pmces="TGLR VGT pmces", tglr_bootstrapped_nonnegative_sampler="TGLR VGT sampler",pmces_uniform_iterative_variant_gene_prior_pip_level_pmces="Iterative VGT pmces", pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped="Iterative VGT sampler")

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=ln_pi_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
    	#scale_fill_manual(values=c("#56B4E9", "darkorchid3"))+
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP=", pip_threshold, " / ", element_class, " / ", new_method))  +
  		geom_hline(yintercept=(1.0-pip_threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12))
  	return(p)
}


make_tgfm_pip_power_plot_varying_eqtl_sample_size_and_prior_method <- function(df, pip_threshold, twas_method, element_class) {

	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(300,500,1000))
	#df$ln_pi_method = factor(df$ln_pi_method, levels=c("uniform", "tglr_variant_gene", "tglr_sparse_variant_gene_tissue", "iterative_variant_gene_tissue", "iterative_variant_gene_tissue_bootstrapped"))
	df$ln_pi_method = factor(df$ln_pi_method, levels=c("uniform", "tglr_bootstrapped_nonnegative_pmces", "tglr_bootstrapped_nonnegative_sampler", "pmces_uniform_iterative_variant_gene_prior_pip_level_pmces", "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"))
	df=df[!is.na(df$ln_pi_method),]	
	df$ln_pi_method = recode(df$ln_pi_method, uniform="Uniform", tglr_bootstrapped_nonnegative_pmces="TGLR VGT pmces", tglr_bootstrapped_nonnegative_sampler="TGLR VGT sampler",pmces_uniform_iterative_variant_gene_prior_pip_level_pmces="Iterative VGT pmces", pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped="Iterative VGT sampler")

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=ln_pi_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("#56B4E9", "darkorchid3"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0("PIP=", pip_threshold, " / ", element_class, " / ", new_method)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12))
  	return(p)
}

make_tgfm_alt_variant_gene_fdr_plot_data_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_nm_variant_alt_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in TGFM gene data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "variant") & (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("gene_only", length(as.character(tmp_df$genetic_element_class))))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("gene", "gene_only", "variant")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))


	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	df2 <- data.frame(method=df$method, eQTL_sample_size=df$eQTL_sample_size, PIP_threshold=rep(pip_threshold, length(df$eQTL_sample_size)), FDR=df$fdr, FDR_lb=df$fdr_lb, FDR_ub=df$fdr_ub)
  	return(df2)

}

make_tgfm_alt_variant_gene_fdr_plot_across_ge_h2 <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_ge_h2s_nm_variant_alt_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GE_h2)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in TGFM gene data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_ge_h2s_gene_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "variant") & (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("gene_only", length(as.character(tmp_df$genetic_element_class))))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GE_h2)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("gene", "gene_only", "variant")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))

	df$eQTL_sample_size = recode(df$eQTL_sample_size, '75'="0.075", '1'="0.1", '5'="0.05")
	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c("0.05", "0.075", "0.1"))


	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[6]
	red_color=brewer.pal(n = 9, name = "Reds")[6]

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0


	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		scale_fill_manual(values=c(red_color, green_color, blue_color))+
  		figure_theme() +
  		labs(x="Gene expression h2", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0-as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")

  	if (plot_expected_fdr) {
  		p = p + geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted')
  	}
  	return(p)

}


make_tgfm_alt_variant_gene_fdr_plot_across_gwas_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_gwas_ss_nm_variant_alt_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GWAS_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in TGFM gene data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_gwas_ss_gene_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "variant") & (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("gene_only", length(as.character(tmp_df$genetic_element_class))))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GWAS_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("gene", "gene_only", "variant")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))



	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[6]
	red_color=brewer.pal(n = 9, name = "Reds")[6]

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0


	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		scale_fill_manual(values=c(red_color, green_color, blue_color))+
  		figure_theme() +
  		labs(x="GWAS sample size", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0-as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")

  	if (plot_expected_fdr) {
  		p = p + geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted')
  	}
  	return(p)

}

make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_nm_variant_alt_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in TGFM gene data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "variant") & (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("gene_only", length(as.character(tmp_df$genetic_element_class))))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("gene", "gene_only", "variant")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))


	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[6]
	red_color=brewer.pal(n = 9, name = "Reds")[6]

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0


	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		scale_fill_manual(values=c(red_color, green_color, blue_color))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0-as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")

  	if (plot_expected_fdr) {
  		p = p + geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted')
  	}
  	return(p)

}


make_tgfm_gene_and_gene_tissue_fdr_plot_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM gene-tissue data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "variant") & (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene-Tissue", length(as.character(tmp_df$genetic_element_class))))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in TGFM gene data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "variant") & (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene", length(as.character(tmp_df$genetic_element_class))))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("Gene-Tissue", "Gene")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)

	#df$method = recode(df$method, variant="Variant", gene="Gene-Tissue")
	df$method = factor(df$method, levels=c("Gene-Tissue", "Gene"))


	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted')  +
  		scale_fill_manual(values=c("dodgerblue3", "palegreen4"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0-as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")
  	return(p)

}


make_tgfm_variant_gene_fdr_plot_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in TGFM gene data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "variant") & (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("gene_only", length(as.character(tmp_df$genetic_element_class))))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("variant", "gene", "gene_only")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))



	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[6]
	red_color=brewer.pal(n = 9, name = "Reds")[6]

	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		scale_fill_manual(values=c(red_color, green_color, blue_color))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0-as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")
  	if (plot_expected_fdr) {
  		p <- p + geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted') 
  	}
  	return(p)

}

make_tgfm_variant_gene_precision_plot_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df$genetic_element_class = as.character(tgfm_calibration_df$genetic_element_class)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$genetic_element_class) != "all") & (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)

	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("variant", "gene")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec)

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=precision, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=precision_lb, ymax=precision_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("#56B4E9", "darkorchid3"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Precision", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")
  	return(p)

}

make_gene_level_not_gene_tissue_fdr_plot_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df = tgfm_calibration_df[as.character(tgfm_calibration_df$genetic_element_class) == "gene",]



	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Load in FOCUS data
	focus_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_pip_", pip_threshold, "_gene_calibration.txt")
	focus_calibration_df <- read.table(focus_calibration_file, header=TRUE)
	focus_calibration_df = focus_calibration_df[as.character(focus_calibration_df$genetic_element_class) == "gene",]
	tmp_df = focus_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in FOCUS data
	focus_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_tg_pip_", pip_threshold, "_gene_calibration.txt")
	focus_calibration_df <- read.table(focus_calibration_file, header=TRUE)
	focus_calibration_df = focus_calibration_df[as.character(focus_calibration_df$genetic_element_class) == "gene",]
	tmp_df = focus_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS-TG", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Load in Coloc data
	coloc_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_coloc_pip_", pip_threshold, "_gene_calibration.txt")
	coloc_calibration_df <- read.table(coloc_calibration_file, header=TRUE)
	coloc_calibration_df = coloc_calibration_df[as.character(coloc_calibration_df$genetic_element_class) == "gene",]
	tmp_df = coloc_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("coloc", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)





	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM", "FOCUS-TG", "FOCUS", "coloc")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)


	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0

	red_color=brewer.pal(n = 9, name = "Reds")[6]
	green_color=brewer.pal(n = 9, name = "Greens")[6]

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		#scale_fill_manual(values=c(brewer.pal(n = 9, name = "Reds")[6], brewer.pal(n = 9, name = "Reds")[5], brewer.pal(n = 9, name = "Reds")[4], brewer.pal(n = 9, name = "Reds")[2]))+
  		scale_fill_manual(values=c(green_color, "#cc5127", "#e4b422", "#ece918"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0-as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")
  	if (plot_expected_fdr) {
  		p = p + geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted')
  	}
  	return(p)

}

make_gene_fdr_plot_data_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df = tgfm_calibration_df[as.character(tgfm_calibration_df$genetic_element_class) == "gene",]


	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Load in FOCUS data
	focus_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_pip_", pip_threshold, "_calibration.txt")
	focus_calibration_df <- read.table(focus_calibration_file, header=TRUE)
	focus_calibration_df = focus_calibration_df[as.character(focus_calibration_df$genetic_element_class) == "gene",]
	tmp_df = focus_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in FOCUS data
	focus_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_tg_pip_", pip_threshold, "_calibration.txt")
	focus_calibration_df <- read.table(focus_calibration_file, header=TRUE)
	focus_calibration_df = focus_calibration_df[as.character(focus_calibration_df$genetic_element_class) == "gene",]
	tmp_df = focus_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS-TG", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Load in Coloc data
	coloc_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_coloc_pip_", pip_threshold, "_calibration.txt")
	coloc_calibration_df <- read.table(coloc_calibration_file, header=TRUE)
	coloc_calibration_df = coloc_calibration_df[as.character(coloc_calibration_df$genetic_element_class) == "gene",]
	tmp_df = coloc_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("coloc", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)





	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM", "FOCUS-TG", "FOCUS", "coloc")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)


	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0


	df2 <- data.frame(method=df$method, eQTL_sample_size=df$eQTL_sample_size,PIP_threshold=rep(pip_threshold, length(df$fdr_ub)), FDR=df$fdr, FDR_lb=df$fdr_lb, FDR_ub=df$fdr_ub)
  	return(df2)
}


make_gene_fdr_plot_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=FALSE) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df = tgfm_calibration_df[as.character(tgfm_calibration_df$genetic_element_class) == "gene",]

	if (FALSE) {
	# Extract data for SuSiE method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_pmces") & (as.character(tgfm_calibration_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("SuSiE", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	}

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	if (FALSE) {
	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM (uniform prior)", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	}

	# Load in FOCUS data
	focus_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_pip_", pip_threshold, "_calibration.txt")
	focus_calibration_df <- read.table(focus_calibration_file, header=TRUE)
	focus_calibration_df = focus_calibration_df[as.character(focus_calibration_df$genetic_element_class) == "gene",]
	tmp_df = focus_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Load in FOCUS data
	focus_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_tg_pip_", pip_threshold, "_calibration.txt")
	focus_calibration_df <- read.table(focus_calibration_file, header=TRUE)
	focus_calibration_df = focus_calibration_df[as.character(focus_calibration_df$genetic_element_class) == "gene",]
	tmp_df = focus_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS-TG", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)


	# Load in Coloc data
	coloc_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_coloc_pip_", pip_threshold, "_calibration.txt")
	coloc_calibration_df <- read.table(coloc_calibration_file, header=TRUE)
	coloc_calibration_df = coloc_calibration_df[as.character(coloc_calibration_df$genetic_element_class) == "gene",]
	tmp_df = coloc_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("coloc", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)





	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM", "FOCUS-TG", "FOCUS", "coloc")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)


	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0

	red_color=brewer.pal(n = 9, name = "Reds")[6]

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		#scale_fill_manual(values=c(brewer.pal(n = 9, name = "Reds")[6], brewer.pal(n = 9, name = "Reds")[5], brewer.pal(n = 9, name = "Reds")[4], brewer.pal(n = 9, name = "Reds")[2]))+
  		scale_fill_manual(values=c(red_color, "#cc5127", "#e4b422", "#ece918"))+
  		#scale_fill_manual(values=c(red_color, "#cc5127", "#ce772a", "#cca22a"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0-as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")
  	if (plot_expected_fdr) {
  		p = p + geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted')
  	}
  	return(p)

}


make_gene_precision_plot_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df = tgfm_calibration_df[as.character(tgfm_calibration_df$genetic_element_class) == "gene",]

	if (FALSE) {
	# Extract data for SuSiE method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_pmces") & (as.character(tgfm_calibration_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("SuSiE", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	}

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)

	if (FALSE) {
	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM (uniform prior)", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	}

	# Load in FOCUS data
	focus_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_pip_", pip_threshold, "_calibration.txt")
	focus_calibration_df <- read.table(focus_calibration_file, header=TRUE)
	focus_calibration_df = focus_calibration_df[as.character(focus_calibration_df$genetic_element_class) == "gene",]
	tmp_df = focus_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)

	# Load in Coloc data
	coloc_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_coloc_pip_", pip_threshold, "_calibration.txt")
	coloc_calibration_df <- read.table(coloc_calibration_file, header=TRUE)
	coloc_calibration_df = coloc_calibration_df[as.character(coloc_calibration_df$genetic_element_class) == "gene",]
	tmp_df = coloc_calibration_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("coloc", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)




	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("coloc", "FOCUS", "TGFM")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec)

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=precision, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=precision_lb, ymax=precision_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("#56B4E9", "darkorchid3"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Precision", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")
  	return(p)

}

make_gene_fdr_plot_across_sample_sizes_for_various_versions_of_tgfm <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	n_detected_vec <- c()
	eQTL_sample_size_vec <- c()
	coverage_vec <- c()
	coverage_lb_vec <- c()
	coverage_ub_vec <- c()
	expected_coverage_vec <- c()

	# Load in TGFM data
	tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	tgfm_calibration_df = tgfm_calibration_df[as.character(tgfm_calibration_df$genetic_element_class) == "gene",]

	
	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_sampler") & (as.character(tgfm_calibration_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM (uniform prior)", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)

	# Extract data for TGFM method
	indices = (as.character(tgfm_calibration_df$twas_method) == "susie_pmces") & (as.character(tgfm_calibration_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_calibration_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM (no sampling, uniform prior)", n_elements))
	n_detected_vec <- c(n_detected_vec, tmp_df$n_detected_elements)
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	coverage_vec <- c(coverage_vec, tmp_df$coverage)
	coverage_lb_vec <- c(coverage_lb_vec, tmp_df$coverage_lb)
	coverage_ub_vec <- c(coverage_ub_vec, tmp_df$coverage_ub)
	expected_coverage_vec <- c(expected_coverage_vec, tmp_df$expected_coverage)



	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM (no sampling, uniform prior)", "TGFM (uniform prior)", "TGFM")), n_detected_elements=n_detected_vec, eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), precision=coverage_vec, precision_ub=coverage_ub_vec, precision_lb=coverage_lb_vec, expected_fdr=1.0-expected_coverage_vec)

	df$fdr = 1.0 - df$precision
	df$fdr_lb = 1.0 - df$precision_ub
	df$fdr_ub = 1.0 - df$precision_lb

	red_color=brewer.pal(n = 9, name = "Reds")[6]
	red_color1=brewer.pal(n = 9, name = "Reds")[2]
	red_color2=brewer.pal(n = 9, name = "Reds")[4]


	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		geom_errorbar(aes(ymin=expected_fdr, ymax=expected_fdr), width=.8, position=position_dodge(.9), linetype='dotted')  +
  		scale_fill_manual(values=c(red_color1, red_color2, red_color))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		geom_hline(yintercept=1.0 - as.numeric(pip_threshold), linetype=2) +
  		theme(legend.position="top")
  	return(p)

}

make_tgfm_gene_and_gene_tissue_power_plot_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM gene tissue data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "variant") & (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene-Tissue", length(as.character(tmp_df$genetic_element_class))))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in TGFM gene data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "variant") & (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene", length(as.character(tmp_df$genetic_element_class))))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)



	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("Gene-Tissue", "Gene")), eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)



	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("dodgerblue3", "palegreen4"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)

}

make_tgfm_variant_gene_power_plot_data_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)

	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)


	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in TGFM gene data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "variant") & (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene_only", length(as.character(tmp_df$genetic_element_class))))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)




	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("variant", "gene", "Gene_only")), eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), PIP_threshold=rep(pip_threshold, length(eQTL_sample_size_vec)), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", Gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))



	return(df)

}

make_tgfm_variant_gene_power_plot_across_ge_h2 <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_ge_h2s_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)

	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)


	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GE_h2)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in TGFM gene data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_ge_h2s_gene_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)

	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "variant") & (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene_only", length(as.character(tmp_df$genetic_element_class))))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GE_h2)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)




	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("variant", "gene", "Gene_only")), eQTL_sample_size=factor(eQTL_sample_size_vec), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", Gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))

	df$eQTL_sample_size = recode(df$eQTL_sample_size, '75'="0.075", '1'="0.1", '5'="0.05")
	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c("0.05", "0.075", "0.1"))


	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[6]
	red_color=brewer.pal(n = 9, name = "Reds")[6]



	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_color, green_color, blue_color))+
  		figure_theme() +
  		labs(x="Gene expression h2", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)

}

make_tgfm_variant_gene_power_plot_across_gwas_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_gwas_ss_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)

	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)


	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GWAS_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in TGFM gene data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_vary_gwas_ss_gene_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)

	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "variant") & (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene_only", length(as.character(tmp_df$genetic_element_class))))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$GWAS_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)




	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("variant", "gene", "Gene_only")), eQTL_sample_size=factor(eQTL_sample_size_vec), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", Gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))

	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[6]
	red_color=brewer.pal(n = 9, name = "Reds")[6]



	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_color, green_color, blue_color))+
  		figure_theme() +
  		labs(x="GWAS sample size", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)

}

make_tgfm_variant_gene_power_plot_across_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)

	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)


	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, as.character(tmp_df$genetic_element_class))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in TGFM gene data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_gene_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df$genetic_element_class = as.character(tgfm_power_df$genetic_element_class)
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$genetic_element_class) != "variant") & (as.character(tgfm_power_df$genetic_element_class) != "all") & (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("Gene_only", length(as.character(tmp_df$genetic_element_class))))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)




	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("variant", "gene", "Gene_only")), eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

	df$method = recode(df$method, variant="TGFM (Variant)", gene="TGFM (Gene-Tissue)", Gene_only="TGFM (Gene)")
	df$method = factor(df$method, levels=c("TGFM (Gene-Tissue)", "TGFM (Gene)", "TGFM (Variant)"))

	blue_color=brewer.pal(n = 9, name = "Blues")[7]
	green_color=brewer.pal(n = 9, name = "Greens")[6]
	red_color=brewer.pal(n = 9, name = "Reds")[6]



	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_color, green_color, blue_color))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)

}

make_gene_power_plot_across_sample_sizes_for_various_versions_of_tgfm <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df = tgfm_power_df[as.character(tgfm_power_df$genetic_element_class) == "gene",]

	
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM (uniform prior)", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$twas_method) == "susie_pmces") & (as.character(tgfm_power_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM (no sampling, uniform prior)", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)



	red_color=brewer.pal(n = 9, name = "Reds")[6]
	red_color1=brewer.pal(n = 9, name = "Reds")[2]
	red_color2=brewer.pal(n = 9, name = "Reds")[4]



	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM (no sampling, uniform prior)", "TGFM (uniform prior)", "TGFM")), eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_color1, red_color2, red_color))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)

}

make_gene_not_gene_tissue_power_plot_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df = tgfm_power_df[as.character(tgfm_power_df$genetic_element_class) == "gene",]


	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)


	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_tg_pip_", pip_threshold, "_gene_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS-TG", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_pip_", pip_threshold, "_gene_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)


	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_coloc_pip_", pip_threshold, "_gene_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("coloc", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)



	red_color=brewer.pal(n = 9, name = "Reds")[6]
	green_color=brewer.pal(n = 9, name = "Greens")[6]

	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM", "FOCUS-TG", "FOCUS", "coloc")), eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(green_color, "#cc5127", "#e4b422", "#ece918"))+
  		#scale_fill_manual(values=c(brewer.pal(n = 9, name = "Reds")[6], brewer.pal(n = 9, name = "Reds")[5], brewer.pal(n = 9, name = "Reds")[4], brewer.pal(n = 9, name = "Reds")[2]))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)

}

make_gene_power_plot_data_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df = tgfm_power_df[as.character(tgfm_power_df$genetic_element_class) == "gene",]

	
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)


	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_tg_pip_", pip_threshold, "_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS-TG", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_pip_", pip_threshold, "_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)


	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_coloc_pip_", pip_threshold, "_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("coloc", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)



	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM", "FOCUS-TG", "FOCUS", "coloc")), eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), PIP_threshold=rep(pip_threshold,length(eQTL_sample_size_vec)), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

  	return(df)

}

make_gene_power_plot_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Initialize vectors for summary df
	method_vec <- c()
	eQTL_sample_size_vec <- c()
	power_vec <- c()
	power_lb_vec <- c()
	power_ub_vec <- c()

	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	tgfm_power_df = tgfm_power_df[as.character(tgfm_power_df$genetic_element_class) == "gene",]

	if (FALSE) {
	# Extract data for SuSiE method
	indices = (as.character(tgfm_power_df$twas_method) == "susie_pmces") & (as.character(tgfm_power_df$ln_pi_method) == "uniform")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("SuSiE", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)
	}
	
	# Extract data for TGFM method
	indices = (as.character(tgfm_power_df$twas_method) == "susie_sampler") & (as.character(tgfm_power_df$ln_pi_method) == "pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped")
	tmp_df = tgfm_power_df[indices,]
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("TGFM", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)


	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_tg_pip_", pip_threshold, "_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS-TG", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)

	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_focus_pip_", pip_threshold, "_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("FOCUS", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)


	# Load in FOCUS data
	focus_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_coloc_pip_", pip_threshold, "_power.txt")
	focus_power_df <- read.table(focus_power_file, header=TRUE)
	focus_power_df = focus_power_df[as.character(focus_power_df$genetic_element_class) == "gene",]
	tmp_df = focus_power_df
	n_elements = dim(tmp_df)[1]
	method_vec <- c(method_vec, rep("coloc", n_elements))
	eQTL_sample_size_vec <- c(eQTL_sample_size_vec, tmp_df$eQTL_sample_size)
	power_vec <- c(power_vec, tmp_df$power)
	power_lb_vec <- c(power_lb_vec, tmp_df$power_lb)
	power_ub_vec <- c(power_ub_vec, tmp_df$power_ub)



	red_color=brewer.pal(n = 9, name = "Reds")[6]
	# Convert into clean data frame
	df <- data.frame(method=factor(method_vec, levels=c("TGFM", "FOCUS-TG", "FOCUS", "coloc")), eQTL_sample_size=factor(eQTL_sample_size_vec, levels=c(300, 500, 1000)), power=power_vec, power_ub=power_ub_vec, power_lb=power_lb_vec)

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c(red_color, "darkgray", "slategray", "#8275ba"))+
  		scale_fill_manual(values=c(red_color, "#cc5127", "#e4b422", "#ece918"))+
  		#scale_fill_manual(values=c(brewer.pal(n = 9, name = "Reds")[6], brewer.pal(n = 9, name = "Reds")[5], brewer.pal(n = 9, name = "Reds")[4], brewer.pal(n = 9, name = "Reds")[2]))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)

}

make_expr_med_fraction_plot_with_standard_errors <- function(df) {
	df$eqtl_sample_size = factor(df$eqtl_sample_size, levels=c(300,500,1000))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=frac_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=frac_h2_lb, ymax=frac_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("#56B4E9", "darkorchid3"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Expected fraction of genetic effects\nmediated through gene expression") +
  		geom_hline(yintercept=.1, linetype=2) 
  	return(p)
}
make_expr_med_fraction_plot_with_standard_errors2 <- function(df) {
	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(300,500,1000))
	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=expression_mediated_fraction)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=expression_mediated_fraction_lb, ymax=expression_mediated_fraction_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("#56B4E9", "darkorchid3"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Expected fraction of TGFM PIP\nmediated through gene expression") +
  		geom_hline(yintercept=.1, linetype=2) 
  	return(p)
}


make_fraction_high_pip_elements_from_gene_expression_plot_with_standard_errors <- function(df) {
	df$eQTL_sample_size = factor(df$eQTL_sample_size, levels=c(300,500,1000))
	df$PIP_threshold = factor(df$PIP_threshold, levels=c(0.1, 0.3, 0.5, 0.7, 0.9))
	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=expression_mediated_fraction, fill=PIP_threshold)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=expression_mediated_fraction_lb, ymax=expression_mediated_fraction_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("#56B4E9", "darkorchid3"))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Expected fraction of fine-mapped genetic\nelements mediated by gene expression", fill="PIP") +
  		geom_hline(yintercept=.1, linetype=2) 
  	return(p)
}


cmp_tgfm_variant_pip_susie_pip_scatter <- function(file_name) {
	df <- read.table(file_name, header=TRUE)

	df$correlated_with_causal_gene = factor(df$correlated_with_causal_gene, levels=c("True", "False"))

	pp <- ggplot(df, aes(y=susie_pip, x=tgfm_pip, color=correlated_with_causal_gene)) +
  		geom_point(size=.1) +
  		figure_theme() +
  		scale_color_manual(values=c("purple", "black")) +
  		labs(y="SuSiE PIP", x="TGFM (Variant) PIP", color="Variant correlated with\nfine-mapped gene") +
  		theme(legend.position="bottom")

  	return(pp)

}


#######################
# Command line args
#######################
global_simulation_name_string = args[1]
simulated_organized_results_dir = args[2]
visualize_simulated_results_dir = args[3]


if (FALSE) {
#####################################################################
#####################################################################
# Genome-wide PRIOR
#####################################################################
#####################################################################
#####################################################################
# Make barplot with standard error showing AVG fraction-mediated per tissue
#####################################################################
version="pmces"
# Load in data
per_tissue_fraction_causal_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_", version, "_uniform_iterative_variant_gene_prior_pip_level_avg_fraction_causal_by_tissue.txt")

per_tissue_fraction_causal_df <- read.table(per_tissue_fraction_causal_file, header=TRUE)
# Make plot
avg_per_tissue_fraction_causal_se_barplot <- make_avg_per_tissue_fraction_causal_se_barplot(per_tissue_fraction_causal_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_iterative_susie_", version, "_avg_fraction_causal_per_tissue.pdf")
ggsave(avg_per_tissue_fraction_causal_se_barplot, file=output_file, width=7.2, height=3.5, units="in")

#####################################################################
# Make barplot with standard error showing Type 1 error across thresholds for bootstrapped iterative VGT PMCES
#####################################################################
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_h2_type_1_error_across_thresholds.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot_across_thresholds(t1e_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_iterative_VGT_PMCES_bootstrapped_med_h2_across_thresholds.pdf")
#ggsave(t1e_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

#####################################################################
# Make barplot with standard error showing Power across thresholds for bootstrapped iterative VGT PMCES
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_mediated_h2_power_across_thresholds.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
power_se_barplot <- make_power_med_h2_se_barplot_across_thresholds(power_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_power_iterative_VGT_PMCES_bootstrapped_med_h2_across_thresholds.pdf")
#ggsave(power_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

# Make joint plot
joint_plot <- plot_grid(t1e_se_barplot + theme(legend.position="none"), power_se_barplot, ncol=1, rel_heights=c(1, 1.4))
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_and_power_iterative_VGT_PMCES_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=6.5, units="in")


#####################################################################
# Make barplot with standard error showing Type 1 error for single threshold for bootstrapped iterative VGT PMCES
#####################################################################
threshold=1e-8
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_h2_type_1_error_across_thresholds.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot_at_single_threshold(t1e_h2_df, threshold)


#####################################################################
# Make barplot with standard error showing Power for single threshold for bootstrapped iterative VGT PMCES
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_mediated_h2_power_across_thresholds.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
power_se_barplot <- make_power_med_h2_se_barplot_at_single_threshold(power_h2_df, threshold)


# Make joint plot
joint_plot <- plot_grid(t1e_se_barplot + theme(legend.position="none"), power_se_barplot+ theme(legend.position="none"), ncol=1, rel_heights=c(1, 1), labels=c("a","b"))
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_and_power_iterative_VGT_PMCES_bootstrapped_med_h2_at_", threshold, ".pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=6.5, units="in")

#####################################################################
# Make barplot with standard error showing Type 1 error using gaussian approximation of bootstrapped iterative VGT PMCES
#####################################################################
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_h2_type_1_error_gaussian_approximation.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
#t1e_se_barplot <- make_type_1_error_med_h2_se_barplot_at_single_threshold(t1e_h2_df, threshold)
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot_gaussian_approximation(t1e_h2_df)


#####################################################################
# Make barplot with standard error showing Power for single threshold for bootstrapped iterative VGT PMCES
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_mediated_h2_power_gaussian_approximation.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
#power_se_barplot <- make_power_med_h2_se_barplot_at_single_threshold(power_h2_df, threshold)
power_se_barplot <- make_power_med_h2_se_barplot_gaussian_approximation(power_h2_df)

# Make joint plot
joint_plot <- plot_grid(t1e_se_barplot + theme(legend.position="none"), power_se_barplot+ theme(legend.position="none"), ncol=1, rel_heights=c(1, 1), labels=c("a","b"))
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_and_power_iterative_VGT_PMCES_bootstrapped_med_h2_gaussian_approximation.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=6.5, units="in")

}



#####################################################################
#####################################################################
# TGFM 
#####################################################################
#####################################################################

if (FALSE) {
#####################################################################
# TGFM Expected fraction of elements mediated by gene expression
#####################################################################
# load in data
expr_med_frac_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_sampler_pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped_expected_fraction_elements_mediated_by_gene_expresssion.txt")
expr_med_frac_df <- read.table(expr_med_frac_file, header=TRUE)
# Make plot
expr_med_frac_se_plot1 <- make_expr_med_fraction_plot_with_standard_errors2(expr_med_frac_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_expected_fraction_TGFM_PIP_mediated_by_expression.pdf")
ggsave(expr_med_frac_se_plot1, file=output_file, width=7.2, height=4.0, units="in")


#####################################################################
# TGFM fraction of high PIP elements mediated by gene expression
#####################################################################
# load in data
expr_med_frac_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_sampler_pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped_fraction_high_pip_elements_mediated_by_gene_expresssion.txt")
expr_med_frac_df <- read.table(expr_med_frac_file, header=TRUE)
# Make plot
expr_med_frac_se_plot2 <- make_fraction_high_pip_elements_from_gene_expression_plot_with_standard_errors(expr_med_frac_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_fraction_TGFM_high_PIP_mediated_by_expression.pdf")
ggsave(expr_med_frac_se_plot2, file=output_file, width=7.2, height=4.0, units="in")


# Joint plot
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_suppfig_estimated_fraction_mediated_by_expression.pdf")
joint_plot <- plot_grid(expr_med_frac_se_plot2, expr_med_frac_se_plot1, ncol=1, labels=c("a","b"))
ggsave(joint_plot, file=output_file, width=7.2, height=7.5, units="in")



#####################################################################
# Make Figure 1
#####################################################################
# Precision plots
pip_threshold <- "0.5"
fdr_plot_5 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
fdr_plot_data_5 <- make_gene_fdr_plot_data_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

pip_threshold <- "0.9"
fdr_plot_9 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
fdr_plot_data_9 <- make_gene_fdr_plot_data_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
power_plot_data_5 <- make_gene_power_plot_data_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

pip_threshold <- "0.9"
power_plot_9 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
power_plot_data_9 <- make_gene_power_plot_data_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

# Extract legend 
legender = get_legend(power_plot_9)

# Make joint plot with cowplot
figure1 <- plot_grid( legender, NULL, plot_grid(fdr_plot_5 +theme(legend.position="none"), fdr_plot_9 +theme(legend.position="none"), power_plot_5 +theme(legend.position="none"), power_plot_9 +theme(legend.position="none"), ncol=2, labels=c("a", "b", "c","d")), ncol=1, rel_heights=c(.05, .03, 1))

# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_figure1.pdf")
ggsave(figure1, file=output_file, width=7.2, height=5.5, units="in")


# Write suppdata for figure 1a-b
fdr_plot_data <- rbind(fdr_plot_data_5, fdr_plot_data_9)
supp_table_file = paste0(visualize_simulated_results_dir, "suppTable_fig1ab.txt")
write.table(fdr_plot_data, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)

# Write suppdata for figure 1c-d
power_plot_data <- rbind(power_plot_data_5, power_plot_data_9)
supp_table_file = paste0(visualize_simulated_results_dir, "suppTable_fig1cd.txt")
write.table(power_plot_data, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)


#####################################################################
# Make Figure 2: TGFM variant and gene precision and power as a function of eQTL sample size
#####################################################################
# Precision plots
pip_threshold <- "0.5"
precision_plot_5 <- make_tgfm_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .6)
pip_threshold <- "0.9"
precision_plot_9 <- make_tgfm_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .6)

# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .25)
pip_threshold <- "0.9"
power_plot_9 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .25)

# Extract legend 
legender = get_legend(power_plot_9)

# Make joint plot with cowplot
joint_figure <- plot_grid( legender, NULL, plot_grid(precision_plot_5 +theme(legend.position="none"), precision_plot_9 +theme(legend.position="none"), power_plot_5 +theme(legend.position="none"), power_plot_9 +theme(legend.position="none"), ncol=2, labels=c("a", "b", "c","d")), ncol=1, rel_heights=c(.05, .03, 1))

# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_figure2_alternative_fdr.pdf")
ggsave(joint_figure, file=output_file, width=7.2, height=5.5, units="in")

# ALT VERSION
# Precision plots
pip_threshold <- "0.5"
precision_plot_5 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .6)
precision_plot_data_5 <- make_tgfm_alt_variant_gene_fdr_plot_data_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
precision_plot_9 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .6)
precision_plot_data_9 <- make_tgfm_alt_variant_gene_fdr_plot_data_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)


# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .25)
power_plot_data_5 <- make_tgfm_variant_gene_power_plot_data_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) 
pip_threshold <- "0.9"
power_plot_9 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0,.25)
power_plot_data_9 <- make_tgfm_variant_gene_power_plot_data_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) 



# Extract legend 
legender = get_legend(power_plot_9)

# Make joint plot with cowplot
joint_figure <- plot_grid( legender, NULL, plot_grid(precision_plot_5 +theme(legend.position="none"), precision_plot_9 +theme(legend.position="none"), power_plot_5 +theme(legend.position="none"), power_plot_9 +theme(legend.position="none"), ncol=2, labels=c("a", "b", "c","d")), ncol=1, rel_heights=c(.05, .03, 1))

# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_figure2.pdf")
ggsave(joint_figure, file=output_file, width=7.2, height=5.5, units="in")


# Write suppdata for figure 2a-b
fdr_plot_data <- rbind(precision_plot_data_5, precision_plot_data_9)
supp_table_file = paste0(visualize_simulated_results_dir, "suppTable_fig2ab.txt")
write.table(fdr_plot_data, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)

# Write suppdata for figure 2c-d
power_plot_data <- rbind(power_plot_data_5, power_plot_data_9)
supp_table_file = paste0(visualize_simulated_results_dir, "suppTable_fig2cd.txt")
write.table(power_plot_data, file=supp_table_file, quote=FALSE, sep="\t", row.names = FALSE)
}




#####################################################################
# Make Supplementary Figure: Versin of Figure 2 where we vary GWAS sample size
#####################################################################
# ALT VERSION
# Precision plots
pip_threshold <- "0.5"
precision_plot_5 <- make_tgfm_alt_variant_gene_fdr_plot_across_gwas_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .6)
pip_threshold <- "0.9"
precision_plot_9 <- make_tgfm_alt_variant_gene_fdr_plot_across_gwas_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .6)


# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_tgfm_variant_gene_power_plot_across_gwas_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .2)
pip_threshold <- "0.9"
power_plot_9 <- make_tgfm_variant_gene_power_plot_across_gwas_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0,.2)



# Extract legend 
legender = get_legend(power_plot_9)

# Make joint plot with cowplot
joint_figure <- plot_grid( legender, NULL, plot_grid(precision_plot_5 +theme(legend.position="none"), precision_plot_9 +theme(legend.position="none"), power_plot_5 +theme(legend.position="none"), power_plot_9 +theme(legend.position="none"), ncol=2, labels=c("a", "b", "c","d")), ncol=1, rel_heights=c(.05, .03, 1))


# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_suppFig_Figure2_vary_gwas_sample_size.pdf")
ggsave(joint_figure, file=output_file, width=7.2, height=5.5, units="in")

#####################################################################
# Make Supplementary Figure: Versin of Figure 2 where we vary GE h2 sample size
#####################################################################
# ALT VERSION
# Precision plots
pip_threshold <- "0.5"
precision_plot_5 <- make_tgfm_alt_variant_gene_fdr_plot_across_ge_h2(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .65)
pip_threshold <- "0.9"
precision_plot_9 <- make_tgfm_alt_variant_gene_fdr_plot_across_ge_h2(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .65)


# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_tgfm_variant_gene_power_plot_across_ge_h2(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .2)
pip_threshold <- "0.9"
power_plot_9 <- make_tgfm_variant_gene_power_plot_across_ge_h2(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0,.2)



# Extract legend 
legender = get_legend(power_plot_9)

# Make joint plot with cowplot
joint_figure <- plot_grid( legender, NULL, plot_grid(precision_plot_5 +theme(legend.position="none"), precision_plot_9 +theme(legend.position="none"), power_plot_5 +theme(legend.position="none"), power_plot_9 +theme(legend.position="none"), ncol=2, labels=c("a", "b", "c","d")), ncol=1, rel_heights=c(.05, .03, 1))


# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_suppFig_Figure2_vary_ge_h2s.pdf")
ggsave(joint_figure, file=output_file, width=7.2, height=5.5, units="in")






if (FALSE) {
#####################################################################
# Make Figure 1-2 hybrid for poster
#####################################################################
# Precision plots
pip_threshold <- "0.5"
fdr_plot_5 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
# Extract legend 
legender = get_legend(power_plot_5 + theme(legend.position="bottom"))
#fig_1_poster <- plot_grid(fdr_plot_5+theme(legend.position="none"), power_plot_5+theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.1))
fig_1_poster <- plot_grid(plot_grid(fdr_plot_5+theme(legend.position="none"), power_plot_5+theme(legend.position="none"), ncol=2), legender, ncol=1, rel_heights=c(1,.13))

pip_threshold <- "0.5"
precision_plot_5 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .6)
# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) + ylim(0, .25)
# Extract legend 
legender = get_legend(power_plot_5+ theme(legend.position="bottom"))
fig_2_poster <- plot_grid(plot_grid(precision_plot_5+theme(legend.position="none"), power_plot_5+theme(legend.position="none"), ncol=2), legender, ncol=1, rel_heights=c(1,.13))


poster_fig <- plot_grid(fig_1_poster, fig_2_poster, ncol=2)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_poster_figure.pdf")
ggsave(poster_fig, file=output_file, width=21.2, height=3.0, units="in")




#####################################################################
# Plot precision over a range of thresholds
#####################################################################
# Precision plots
pip_threshold <- "0.3"
fdr_plot_3 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.5"
fdr_plot_5 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.7"
fdr_plot_7 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.9"
fdr_plot_9 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.95"
fdr_plot_95 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.99"
fdr_plot_99 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)


# Extract legend 
legender = get_legend(fdr_plot_99)

# Make joint plot with cowplot
figure <- plot_grid(legender, plot_grid(fdr_plot_3+theme(legend.position="none"), fdr_plot_5+theme(legend.position="none"), fdr_plot_7+theme(legend.position="none"), fdr_plot_9+theme(legend.position="none"), fdr_plot_95+theme(legend.position="none"), fdr_plot_99+theme(legend.position="none"),ncol=2, labels=c("a", "b", "c","d", "e", "f")),ncol=1, rel_heights=c(.05,1))


# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_method_precision_pip_range.pdf")
ggsave(figure, file=output_file, width=7.2, height=5.5, units="in")

# Precision plots
pip_threshold <- "0.3"
fdr_plot_3 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.5"
fdr_plot_5 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.7"
fdr_plot_7 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.9"
fdr_plot_9 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.95"
fdr_plot_95 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)
pip_threshold <- "0.99"
fdr_plot_99 <- make_tgfm_alt_variant_gene_fdr_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold,plot_expected_fdr=TRUE)


# Extract legend 
legender = get_legend(fdr_plot_99)

# Make joint plot with cowplot
figure <- plot_grid(legender, plot_grid(fdr_plot_3+theme(legend.position="none"), fdr_plot_5+theme(legend.position="none"), fdr_plot_7+theme(legend.position="none"), fdr_plot_9+theme(legend.position="none"), fdr_plot_95+theme(legend.position="none"), fdr_plot_99+theme(legend.position="none"),ncol=2, labels=c("a", "b", "c","d", "e", "f")),ncol=1, rel_heights=c(.05,1))


# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_tissue_gene_variant_precision_pip_range.pdf")
ggsave(figure, file=output_file, width=7.2, height=5.5, units="in")




#####################################################################
# Plot power over a range of thresholds
#####################################################################
# Precision plots
pip_threshold <- "0.3"
precision_plot_3 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.5"
precision_plot_5 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.7"
precision_plot_7 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
precision_plot_9 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.95"
precision_plot_95 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.99"
precision_plot_99 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)


# Extract legend 
legender = get_legend(precision_plot_99)

# Make joint plot with cowplot
figure <- plot_grid(legender, plot_grid(precision_plot_3+theme(legend.position="none"), precision_plot_5+theme(legend.position="none"), precision_plot_7+theme(legend.position="none"), precision_plot_9+theme(legend.position="none"), precision_plot_95+theme(legend.position="none"), precision_plot_99+theme(legend.position="none"),ncol=2, labels=c("a", "b", "c","d", "e", "f")),ncol=1, rel_heights=c(.05,1))

# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_method_power_pip_range.pdf")
ggsave(figure, file=output_file, width=7.2, height=5.5, units="in")

# Precision plots
pip_threshold <- "0.3"
precision_plot_3 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.5"
precision_plot_5 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.7"
precision_plot_7 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
precision_plot_9 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.95"
precision_plot_95 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.99"
precision_plot_99 <- make_tgfm_variant_gene_power_plot_across_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)


# Extract legend 
legender = get_legend(precision_plot_99)

# Make joint plot with cowplot
figure <- plot_grid(legender, plot_grid(precision_plot_3+theme(legend.position="none"), precision_plot_5+theme(legend.position="none"), precision_plot_7+theme(legend.position="none"), precision_plot_9+theme(legend.position="none"), precision_plot_95+theme(legend.position="none"), precision_plot_99+theme(legend.position="none"),ncol=2, labels=c("a", "b", "c","d", "e", "f")),ncol=1, rel_heights=c(.05,1))


# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_tissue_gene_variant_power_pip_range.pdf")
ggsave(figure, file=output_file, width=7.2, height=5.5, units="in")


#####################################################################
# Make version of Figure 1 comparing TGFM with and with out prior
#####################################################################
# Precision plots
pip_threshold <- "0.5"
fdr_plot_5 <- make_gene_fdr_plot_across_sample_sizes_for_various_versions_of_tgfm(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
fdr_plot_9 <- make_gene_fdr_plot_across_sample_sizes_for_various_versions_of_tgfm(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_gene_power_plot_across_sample_sizes_for_various_versions_of_tgfm(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
power_plot_9 <- make_gene_power_plot_across_sample_sizes_for_various_versions_of_tgfm(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

# Extract legend 
legender = get_legend(power_plot_9)

# Make joint plot with cowplot
figure <- plot_grid( legender, NULL, plot_grid(fdr_plot_5 +theme(legend.position="none"), fdr_plot_9 +theme(legend.position="none"), power_plot_5 +theme(legend.position="none"), power_plot_9 +theme(legend.position="none"), ncol=2, labels=c("a", "b", "c","d")), ncol=1, rel_heights=c(.05, .03, 1))

# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_various_tgfm_models_precision_and_power.pdf")
ggsave(figure, file=output_file, width=7.2, height=5.5, units="in")




#####################################################################
# Make gene-level (not gene-tissue) version of Figure 1
#####################################################################
# Precision plots
pip_threshold <- "0.5"
fdr_plot_5 <- make_gene_level_not_gene_tissue_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=TRUE)
pip_threshold <- "0.9"
fdr_plot_9 <- make_gene_level_not_gene_tissue_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, plot_expected_fdr=TRUE)

# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_gene_not_gene_tissue_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
power_plot_9 <- make_gene_not_gene_tissue_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

# Extract legend 
legender = get_legend(power_plot_9)

# Make joint plot with cowplot
figure1 <- plot_grid( legender, NULL, plot_grid(fdr_plot_5 +theme(legend.position="none"), fdr_plot_9 +theme(legend.position="none"), power_plot_5 +theme(legend.position="none"), power_plot_9 +theme(legend.position="none"), ncol=2, labels=c("a", "b", "c","d")), ncol=1, rel_heights=c(.05, .03, 1))

# Make joint plot
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_level_version_of_figure1.pdf")
ggsave(figure1, file=output_file, width=7.2, height=5.5, units="in")
}



if (FALSE) {
#####################################################################
# Scatter plot comparing TGFM (Variant) PIP with SuSiE PIP
#####################################################################
tgfm_variant_susie_cmp_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_tgfm_variant_pip_susie_variant_pip_comparison.txt")
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_variant_susie_variant_pip_cmp_scatter.pdf")
scatter <- cmp_tgfm_variant_pip_susie_pip_scatter(tgfm_variant_susie_cmp_file)
ggsave(scatter, file=output_file, width=7.2, height=4.6, units="in")
}











































# REcent old
if (FALSE) {
#####################################################################
# Make barplot with standard error showing AVG fraction-mediated by gene expression according to prior
#####################################################################
# load in data
expr_med_frac_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_avg_expected_fraction_expression_mediated_disease_components.txt")
expr_med_frac_df <- read.table(expr_med_frac_file, header=TRUE)
# Make plot
expr_med_frac_se_plot <- make_expr_med_fraction_plot_with_standard_errors(expr_med_frac_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_iterative_VGT_PMCES_bootstrapped_prior_expectation_of_expr_mediated_fraction_across_sample_sizes.pdf")
ggsave(expr_med_frac_se_plot, file=output_file, width=7.2, height=4.0, units="in")

}




















































if (FALSE) {
#####################################################################
# Genome-wide analysis
#####################################################################
#####################################################################
# Make barplot with standard error showing AVG heritability estimates
#####################################################################
# load in data
avg_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_avg_total_est_h2.txt")
avg_h2_df <- read.table(avg_h2_file, header=TRUE)
# Make plot
avg_h2_se_barplot <- make_avg_h2_se_barplot(avg_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_avg_total_h2_se_barplot.pdf")
ggsave(avg_h2_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

#####################################################################
# Make barplot with standard error showing AVG fraction of heritability mediated by gene expression in any tissue
#####################################################################
# load in data
avg_frac_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_avg_fraction_expression_med_h2.txt")
avg_frac_h2_df <- read.table(avg_frac_h2_file, header=TRUE)
# Make plot
avg_frac_h2_se_barplot <- make_avg_fraction_h2_se_barplot(avg_frac_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_avg_fraction_h2_mediated_se_barplot.pdf")
ggsave(avg_frac_h2_se_barplot, file=output_file, width=7.2, height=4.5, units="in")


#####################################################################
# Make barplot with standard error showing AVG heritability estimates per tissue
#####################################################################
# Load in data
per_tissue_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_avg_mediated_h2_by_tissue.txt")
per_tissue_h2_df <- read.table(per_tissue_h2_file, header=TRUE)
# Make plot
avg_per_tissue_h2_se_barplot <- make_avg_per_tissue_h2_se_barplot(per_tissue_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_avg_h2_per_tissue.pdf")
ggsave(avg_per_tissue_h2_se_barplot, file=output_file, width=7.2, height=4.5, units="in")


#####################################################################
# Make barplot with standard error showing AVG heritability estimates per tissue given sparse estimates
#####################################################################
# Load in data
per_tissue_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_avg_mediated_h2_by_tissue_sparse_est.txt")
per_tissue_h2_df <- read.table(per_tissue_h2_file, header=TRUE)
# Make plot
avg_per_tissue_h2_se_barplot <- make_avg_per_tissue_h2_se_barplot(per_tissue_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_avg_h2_per_tissue_sparse_est.pdf")
ggsave(avg_per_tissue_h2_se_barplot, file=output_file, width=7.2, height=4.5, units="in")


#####################################################################
# Make barplot with standard error showing AVG tau estimates per tissue given nonnegative estimates
#####################################################################
# Load in data
per_tissue_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_full_anno_avg_mediated_nonnegative_tau_tissue_est.txt")
per_tissue_h2_df <- read.table(per_tissue_h2_file, header=TRUE)
# Make plot
avg_per_tissue_h2_se_barplot <- make_avg_per_tissue_tau_se_barplot(per_tissue_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_avg_tau_per_tissue_nonnegative_tglr_est.pdf")
ggsave(avg_per_tissue_h2_se_barplot, file=output_file, width=7.2, height=3.5, units="in")
}


if (FALSE) {
version="sampler"
# Load in data
per_tissue_fraction_causal_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_", version, "_uniform_iterative_variant_gene_prior_pip_level_avg_fraction_causal_by_tissue.txt")

per_tissue_fraction_causal_df <- read.table(per_tissue_fraction_causal_file, header=TRUE)
# Make plot
avg_per_tissue_fraction_causal_se_barplot <- make_avg_per_tissue_fraction_causal_se_barplot(per_tissue_fraction_causal_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_iterative_susie_", version, "_avg_fraction_causal_per_tissue.pdf")
ggsave(avg_per_tissue_fraction_causal_se_barplot, file=output_file, width=7.2, height=3.5, units="in")
}


if (FALSE) {
#####################################################################
# Make barplot with standard error showing Power
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_full_anno_mediated_h2_power.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
power_se_barplot <- make_power_med_h2_se_barplot(power_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_power_med_h2.pdf")
ggsave(power_se_barplot, file=output_file, width=7.2, height=4.5, units="in")


#####################################################################
# Make barplot with standard error showing Type 1 error
#####################################################################
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_full_anno_mediated_h2_type_1_error.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot(t1e_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_med_h2.pdf")
ggsave(t1e_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

#####################################################################
# Make barplot with standard error showing Power with sparse model
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_sparse_mediated_h2_power.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
power_se_barplot <- make_power_med_h2_se_barplot(power_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_power_sparse_med_h2.pdf")
ggsave(power_se_barplot, file=output_file, width=7.2, height=4.5, units="in")


#####################################################################
# Make barplot with standard error showing Type 1 error with sparse model
#####################################################################
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_sparse_mediated_h2_type_1_error.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot(t1e_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_sparse_med_h2.pdf")
ggsave(t1e_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

#####################################################################
# Make barplot with standard error showing Type 1 error across thresholds for bootstrapped non-negative TGLR
#####################################################################
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_full_anno_nonnegative_bootstrapped_mediated_h2_type_1_error_across_thresholds.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot_across_thresholds(t1e_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_tglr_nonnegative_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(t1e_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

#####################################################################
# Make barplot with standard error showing Power across thresholds for bootstrapped non-negative TGLR
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_full_anno_nonnegative_bootstrapped_mediated_h2_power_across_thresholds.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
power_se_barplot <- make_power_med_h2_se_barplot_across_thresholds(power_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_power_tglr_nonnegative_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(power_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

# Make joint plot
joint_plot <- plot_grid(t1e_se_barplot, power_se_barplot, ncol=1)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_and_power_tglr_nonnegative_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=7.5, units="in")


#####################################################################
# Make barplot with standard error showing Type 1 error across thresholds for bootstrapped non-negative TGLR (genotype intercept)
#####################################################################
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_genotype_intercept_nonnegative_bootstrapped_mediated_h2_type_1_error_across_thresholds.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot_across_thresholds(t1e_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_tglr_genotype_intercept_nonnegative_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(t1e_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

#####################################################################
# Make barplot with standard error showing Power across thresholds for bootstrapped non-negative TGLR (genotype intercept)
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_pmces_genotype_intercept_nonnegative_bootstrapped_mediated_h2_power_across_thresholds.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
power_se_barplot <- make_power_med_h2_se_barplot_across_thresholds(power_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_power_tglr_genotype_intercept_nonnegative_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(power_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

# Make joint plot
joint_plot <- plot_grid(t1e_se_barplot, power_se_barplot, ncol=1)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_and_power_tglr_genotype_intercept_nonnegative_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=7.5, units="in")




#####################################################################
# Make barplot with standard error showing Type 1 error across thresholds for bootstrapped iterative VGT Sampler
#####################################################################
# load in data
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_sampler_uniform_iterative_variant_gene_prior_pip_level_h2_type_1_error_across_thresholds.txt")
t1e_h2_df <- read.table(t1e_h2_file, header=TRUE)
# Make plot
t1e_se_barplot <- make_type_1_error_med_h2_se_barplot_across_thresholds(t1e_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_iterative_VGT_sampler_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(t1e_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

#####################################################################
# Make barplot with standard error showing Power across thresholds for bootstrapped iterative VGT Sampler
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_susie_sampler_uniform_iterative_variant_gene_prior_pip_level_mediated_h2_power_across_thresholds.txt")
power_h2_df <- read.table(power_h2_file, header=TRUE)
# Make plot
power_se_barplot <- make_power_med_h2_se_barplot_across_thresholds(power_h2_df)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_power_iterative_VGT_sampler_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(power_se_barplot, file=output_file, width=7.2, height=4.5, units="in")

# Make joint plot
joint_plot <- plot_grid(t1e_se_barplot, power_se_barplot, ncol=1)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_type_1_error_and_power_iterative_VGT_sampler_bootstrapped_med_h2_across_thresholds.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=7.5, units="in")

}



if (FALSE) {
#####################################################################
# TGFM
#####################################################################
## Quantities of interest
## 1. variant vs gene
## 2. eQTL sample size
## 3. PMCES vs Sampler
## 4. Prior choice

#####################################################################
# Calibration at varying pip thresholds (plots) stratefied by sample size (x-axis) and variant vs gene
#####################################################################
# PMCES
ln_pi_method="uniform"
twas_method="susie_pmces"

pip_threshold_arr <- c("0.5", "0.7", "0.9")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
	calibration_df = calibration_df[as.character(calibration_df$twas_method) == twas_method,]
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) != "all",]

	if (twas_method == "susie_pmces") {
		new_method = "TGFM_pmces"
	}
	if (twas_method == "susie_sampler") {
		new_method = "TGFM_sampler"
	}

	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_element_class(calibration_df, as.numeric(pip_threshold), ln_pi_method, new_method)
	plots[[pip_iter]] = calibration_barplot

	if (pip_threshold=="0.9") {
		output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", twas_method, "_", ln_pi_method, "_fdr_vary_eqtl_sample_size_and_element_class_pip_", pip_threshold, ".pdf")
		ggsave(calibration_barplot, file=output_file, width=7.2, height=3.7, units="in")
	}

}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", twas_method, "_", ln_pi_method, "_fdr_vary_eqtl_sample_size_and_element_class.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")

# SAMPLER
ln_pi_method="uniform"
twas_method="susie_sampler"

pip_threshold_arr <- c("0.5", "0.7", "0.9")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
	calibration_df = calibration_df[as.character(calibration_df$twas_method) == twas_method,]
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) != "all",]

	if (twas_method == "susie_pmces") {
		new_method = "TGFM_pmces"
	}
	if (twas_method == "susie_sampler") {
		new_method = "TGFM_sampler"
	}

	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_element_class(calibration_df, as.numeric(pip_threshold), ln_pi_method, new_method)
	plots[[pip_iter]] = calibration_barplot
	if (pip_threshold=="0.9") {
		output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", twas_method, "_", ln_pi_method, "_fdr_vary_eqtl_sample_size_and_element_class_pip_", pip_threshold, ".pdf")
		ggsave(calibration_barplot, file=output_file, width=7.2, height=3.7, units="in")
	}
}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", twas_method, "_", ln_pi_method, "_fdr_vary_eqtl_sample_size_and_element_class.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")




#####################################################################
# Calibration at varying pip thresholds (plots) while just showing gene FDR stratefied by sample size (x-axis) and pmces vs sampler
#####################################################################
# Uniform Prior
ln_pi_method="uniform"
element_class="gene"
pip_threshold_arr <- c("0.5", "0.7", "0.9")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) == element_class,]

	calibration_df$twas_method = recode(calibration_df$twas_method, susie_pmces="TGFM_pmces", susie_sampler="TGFM_sampler")

	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_twas_method(calibration_df, as.numeric(pip_threshold), ln_pi_method, element_class)
	plots[[pip_iter]] = calibration_barplot
	if (pip_threshold=="0.9") {
		output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", ln_pi_method, "_prior_fdr_vary_eqtl_sample_size_and_twas_method_pip_", pip_threshold, ".pdf")
		ggsave(calibration_barplot, file=output_file, width=7.2, height=3.7, units="in")
	}

}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", ln_pi_method, "_prior_fdr_vary_eqtl_sample_size_and_twas_method.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")

# TGLR variant-gene prior
ln_pi_method="tglr_variant_gene"
element_class="gene"
pip_threshold_arr <- c("0.5", "0.7", "0.9")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$ln_pi_method) == ln_pi_method,]
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) == element_class,]

	calibration_df$twas_method = recode(calibration_df$twas_method, susie_pmces="TGFM_pmces", susie_sampler="TGFM_sampler")

	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_twas_method(calibration_df, as.numeric(pip_threshold), ln_pi_method, element_class)
	plots[[pip_iter]] = calibration_barplot
}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", ln_pi_method, "_prior_fdr_vary_eqtl_sample_size_and_twas_method.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")


#####################################################################
# Calibration at varying pip thresholds (plots) while just showing gene FDR stratefied by sample size (x-axis) and prior methods
#####################################################################
# PMCES
twas_method="susie_pmces"
element_class="gene"
pip_threshold_arr <- c("0.5", "0.7", "0.9", "0.99")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) == element_class,]
	calibration_df = calibration_df[as.character(calibration_df$twas_method) == twas_method,]

	if (twas_method == "susie_pmces") {
		new_method = "TGFM_pmces"
	}
	if (twas_method == "susie_sampler") {
		new_method = "TGFM_sampler"
	}


	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_prior_method(calibration_df, as.numeric(pip_threshold), new_method, element_class)
	plots[[pip_iter]] = calibration_barplot
}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", twas_method, "_fdr_vary_eqtl_sample_size_and_prior_method.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")

# Sampler
twas_method="susie_sampler"
element_class="gene"
pip_threshold_arr <- c("0.5", "0.7", "0.9", "0.99")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
	calibration_df <- read.table(calibration_file, header=TRUE)
	calibration_df = calibration_df[as.character(calibration_df$genetic_element_class) == element_class,]
	calibration_df = calibration_df[as.character(calibration_df$twas_method) == twas_method,]

	if (twas_method == "susie_pmces") {
		new_method = "TGFM_pmces"
	}
	if (twas_method == "susie_sampler") {
		new_method = "TGFM_sampler"
	}


	# Make plot
	calibration_barplot <- make_tgfm_pip_fdr_plot_varying_eqtl_sample_size_and_prior_method(calibration_df, as.numeric(pip_threshold), new_method, element_class)
	plots[[pip_iter]] = calibration_barplot
	if (pip_threshold=="0.9") {
		output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", twas_method, "_fdr_vary_eqtl_sample_size_and_prior_method_pip_", pip_threshold, ".pdf")
		ggsave(calibration_barplot, file=output_file, width=7.2, height=3.7, units="in")
	}

}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", twas_method, "_fdr_vary_eqtl_sample_size_and_prior_method.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")







#####################################################################
# Power at varying pip thresholds (plots) while just showing gene FDR stratefied by sample size (x-axis) and prior methods
#####################################################################
# Gene
twas_method="susie_sampler"
element_class="gene"
pip_threshold_arr <- c("0.5", "0.7", "0.9", "0.99")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	power_df <- read.table(power_file, header=TRUE)
	power_df = power_df[as.character(power_df$genetic_element_class) == element_class,]
	power_df = power_df[as.character(power_df$twas_method) == twas_method,]

	if (twas_method == "susie_pmces") {
		new_method = "TGFM_pmces"
	}
	if (twas_method == "susie_sampler") {
		new_method = "TGFM_sampler"
	}

	power_barplot <- make_tgfm_pip_power_plot_varying_eqtl_sample_size_and_prior_method(power_df, pip_threshold, new_method, element_class)
	plots[[pip_iter]] = power_barplot
	if (pip_threshold=="0.9") {
		output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", twas_method, "_power_vary_eqtl_sample_size_and_prior_method_pip_", pip_threshold,".pdf")
		ggsave(power_barplot, file=output_file, width=7.2, height=3.7, units="in")
	}
}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", twas_method, "_power_vary_eqtl_sample_size_and_prior_method.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")

# Variant
twas_method="susie_sampler"
element_class="variant"
pip_threshold_arr <- c("0.5", "0.7", "0.9", "0.99")
plots <- list()
for (pip_iter in 1:length(pip_threshold_arr)) {
	pip_threshold=pip_threshold_arr[pip_iter]
	# Load in data
	power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_power.txt")
	power_df <- read.table(power_file, header=TRUE)
	power_df = power_df[as.character(power_df$genetic_element_class) == element_class,]
	power_df = power_df[as.character(power_df$twas_method) == twas_method,]

	if (twas_method == "susie_pmces") {
		new_method = "TGFM_pmces"
	}
	if (twas_method == "susie_sampler") {
		new_method = "TGFM_sampler"
	}

	power_barplot <- make_tgfm_pip_power_plot_varying_eqtl_sample_size_and_prior_method(power_df, pip_threshold, new_method, element_class)
	plots[[pip_iter]] = power_barplot
}
joint <- plot_grid(plotlist=plots, ncol=1)
# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_", element_class, "_", twas_method, "_power_vary_eqtl_sample_size_and_prior_method.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")


}

