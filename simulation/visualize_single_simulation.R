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



#######################
# Command line args
#######################
global_simulation_name_string = args[1]
simulated_organized_results_dir = args[2]
visualize_simulated_results_dir = args[3]

if (FALSE) {
#####################################################################
# Fraction of genes detected
#####################################################################
# Load in data
number_detected_genes_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_fraction_of_detected_heritable_genes.txt")
n_genes_df <- read.table(number_detected_genes_file, header=TRUE)
# Output file
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_number_of_susie_detected_genes.pdf")
# Make plot
n_genes_se_barplot <- make_n_detected_genes_se_barplot(n_genes_df)
ggsave(n_genes_se_barplot, file=output_file, width=7.2, height=3.2, units="in")

#####################################################################
# TGFM-SLDSC
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
# Make barplot with standard error showing Power
#####################################################################
# load in data
power_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_mediated_h2_power.txt")
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
t1e_h2_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_mediated_h2_type_1_error.txt")
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

}




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
eqtl_sample_size = eqtl_ss_arr[eqtl_ss_iter]
# Load in data
calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string,"_tgfm_pip_", pip_threshold, "_calibration.txt")
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
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_fdr_vary_twas_method_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")

#####################################################################
# Power at PIP threshold of .9 at twas versions for variant and gene assuming single prior and fixed sample size
#####################################################################

pip_threshold=0.9
ln_pi_method="uniform"
initialization_version="best"
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



# Make plot
power_barplot <- make_tgfm_pip_power_plot_varying_twas_method_and_element_class(power_df, pip_threshold, ln_pi_method, eqtl_sample_size)
plots[[eqtl_ss_iter]] = power_barplot


}
joint <- plot_grid(plotlist=plots, ncol=2)

# Save to output
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_tgfm_power_vary_twas_method_with_", ln_pi_method, "_prior_", pip_threshold, "_pip.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.7, units="in")



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







