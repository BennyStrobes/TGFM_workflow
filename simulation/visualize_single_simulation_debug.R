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


make_gene_fdr_plot_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, gene_level=FALSE, plot_expected_fdr=FALSE) {

	# Load in TGFM data
	if (gene_level == FALSE) {
		tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_gt_arch_1_caus_t","_tgfm_pip_", pip_threshold, "_missing_causal_tissue_calibration.txt")
	} else {
		tgfm_calibration_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_gt_arch_1_caus_t","_tgfm_pip_", pip_threshold, "_missing_causal_tissue_gene_calibration.txt")
	}

	tgfm_calibration_df <- read.table(tgfm_calibration_file, header=TRUE)
	df = tgfm_calibration_df[as.character(tgfm_calibration_df$genetic_element_class) == "gene",]




	# Convert into clean data frame


	df$fdr = 1.0 - df$coverage
	df$fdr_lb = 1.0 - df$coverage_ub
	df$fdr_ub = 1.0 - df$coverage_lb

	df$fdr_lb[df$fdr_lb < 0.0] = 0.0

	df$eQTL_sample_size = factor(df$eQTL_sample_size)

	#red_color=brewer.pal(n = 9, name = "Reds")[6]

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=fdr, fill=missingness_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.3, position=position_dodge(.9), size=.5)  +
  		#scale_fill_manual(values=c(brewer.pal(n = 9, name = "Reds")[6], brewer.pal(n = 9, name = "Reds")[5], brewer.pal(n = 9, name = "Reds")[4], brewer.pal(n = 9, name = "Reds")[2]))+
  		#scale_fill_manual(values=c(red_color, "#cc5127", "#ece918"))+
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


make_gene_power_plot_across_methods_and_sample_sizes <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	# Load in TGFM data
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_gt_arch_1_caus_t","_tgfm_pip_", pip_threshold, "_missing_causal_tissue_power.txt")
	tgfm_power_df <- read.table(tgfm_power_file, header=TRUE)
	df = tgfm_power_df[as.character(tgfm_power_df$genetic_element_class) == "gene",]

	df$eQTL_sample_size = factor(df$eQTL_sample_size)

	# Convert into clean data frame

	p<-ggplot(data=df, aes(x=eQTL_sample_size, y=power, fill=missingness_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c(red_color, "darkgray", "slategray", "#8275ba"))+
  		#scale_fill_manual(values=c(red_color, "#cc5127", "#e4b422", "#ece918"))+
  		#scale_fill_manual(values=c(brewer.pal(n = 9, name = "Reds")[6], brewer.pal(n = 9, name = "Reds")[5], brewer.pal(n = 9, name = "Reds")[4], brewer.pal(n = 9, name = "Reds")[2]))+
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="", title=paste0("PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)
}

make_stratefied_power_plot_single_ss <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss) {
	tgfm_power_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_gt_arch_1_caus_t","_tgfm_pip_", pip_threshold, "_missing_causal_tissue_causal_gene_tissue_pair_stratification.txt")

	df<- read.table(tgfm_power_file, header=TRUE,sep="\t")

	df = df[as.character(df$eQTL_sample_size)==eqtl_ss,]

	df$genetic_element_class = factor(df$genetic_element_class, levels=c("causal_gt", "best_tagging_gt", "other_gt", "nm_variant"))

	p<-ggplot(data=df, aes(x=genetic_element_class, y=power, fill=missingness_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c(red_color, "darkgray", "slategray", "#8275ba"))+
  		#scale_fill_manual(values=c(red_color, "#cc5127", "#e4b422", "#ece918"))+
  		#scale_fill_manual(values=c(brewer.pal(n = 9, name = "Reds")[6], brewer.pal(n = 9, name = "Reds")[5], brewer.pal(n = 9, name = "Reds")[4], brewer.pal(n = 9, name = "Reds")[2]))+
  		figure_theme() +
  		labs(x="", y="Power", fill="", title=paste0("eQTL SS = ", eqtl_ss, " / PIP >= ", pip_threshold)) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) + 
  		theme(legend.position="top")
  	return(p)


}

make_stratefied_power_plot <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold) {
	eqtl_ss <- "300"
	p1 <- make_stratefied_power_plot_single_ss(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss)
	eqtl_ss <- "500"
	p2 <- make_stratefied_power_plot_single_ss(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss)
	eqtl_ss <- "1000"
	p3 <- make_stratefied_power_plot_single_ss(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss)

	legender <- get_legend(p1)
	joint_plot <- plot_grid(p1 + theme(legend.position="none"), p2+ theme(legend.position="none"), p3+ theme(legend.position="none"), ncol=1)

	joint_plot2 <- plot_grid(joint_plot, legender, rel_heights=c(1,.1), ncol=1)
	return(joint_plot2)
}


make_correlation_with_causal_gene_scatter <- function(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss) {
	tgfm_corr_file <- paste0(simulated_organized_results_dir, "organized_simulation_", global_simulation_name_string, "_gt_arch_1_caus_t","_tgfm_pip_", pip_threshold,"_not_best_tagging_gene_tissue_pairs_when_missing_causal_tissue.txt")
	df <- read.table(tgfm_corr_file, header=TRUE, sep="\t")
	df = df[as.character(df$eQTL_sample_size) == eqtl_ss,]

	p <- ggplot(df, aes(y=correlation_with_causal_gene, x=best_correlation_with_causal_gene, color=PIP)) +
  		geom_point() +
  		figure_theme() +
  		geom_abline() +
  		labs(title=paste0("eQTL SS: ", eqtl_ss))

  	return(p)



}




#######################
# Command line args
#######################
global_simulation_name_string = args[1]
simulated_organized_results_dir = args[2]
visualize_simulated_results_dir = args[3]



####### Scatter plot showing diff between top curr and best corr
eqtl_ss <- "300"
pip_threshold <- "0.5"
p1 <- make_correlation_with_causal_gene_scatter(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss)

eqtl_ss <- "500"
pip_threshold <- "0.5"
p2 <- make_correlation_with_causal_gene_scatter(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss)

eqtl_ss <- "1000"
pip_threshold <- "0.5"
p3 <- make_correlation_with_causal_gene_scatter(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, eqtl_ss)

joint_plot <- plot_grid(p1,p2,p3,ncol=1)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_corr_with_causal_gene_scatter.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=7.5, units="in")



if (FALSE) {
# FDR plots
pip_threshold <- "0.5"
fdr_plot_5 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
fdr_plot_9 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)


# Power plots
pip_threshold <- "0.5"
power_plot_5 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)
pip_threshold <- "0.9"
power_plot_9 <- make_gene_power_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

legender <- get_legend(power_plot_9)

joint_plot_fdr <- plot_grid(fdr_plot_5 + theme(legend.position="none"), fdr_plot_9+ theme(legend.position="none"), ncol=2)
joint_plot_power <- plot_grid(power_plot_5 + theme(legend.position="none"), power_plot_9+ theme(legend.position="none"), ncol=2)
joint_plot <- plot_grid(joint_plot_fdr, joint_plot_power, legender, rel_heights=c(1,1,.2), ncol=1)
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_missing_causal_tissue_fdr_and_power.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=6.5, units="in")


# Gene calibration
# FDR plots
pip_threshold <- "0.5"
fdr_plot_5 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, gene_level=TRUE)
pip_threshold <- "0.9"
fdr_plot_9 <- make_gene_fdr_plot_across_methods_and_sample_sizes(simulated_organized_results_dir, global_simulation_name_string, pip_threshold, gene_level=TRUE)
legender <- get_legend(fdr_plot_5)
joint_plot_fdr <- plot_grid(fdr_plot_5 + theme(legend.position="none"), fdr_plot_9+ theme(legend.position="none"), ncol=2)
joint_plot <- plot_grid(joint_plot_fdr, legender, ncol=1, rel_heights=c(1,.2))
output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_missing_causal_tissue_gene_level_fdr.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=4.2, units="in")


# Causal gene-tissue stratefied power
pip_threshold <- "0.5"
stratefied_power_5 <- make_stratefied_power_plot(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_tissue_stratefied_power_5.pdf")
ggsave(stratefied_power_5, file=output_file, width=7.2, height=6.2, units="in")

pip_threshold <- "0.9"
stratefied_power_9 <- make_stratefied_power_plot(simulated_organized_results_dir, global_simulation_name_string, pip_threshold)

output_file <- paste0(visualize_simulated_results_dir, "simulation_", global_simulation_name_string, "_gene_tissue_stratefied_power_9.pdf")
ggsave(stratefied_power_9, file=output_file, width=7.2, height=6.2, units="in")
}