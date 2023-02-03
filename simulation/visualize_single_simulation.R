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
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=power)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Power") +
  		theme(legend.position="bottom") 
  	return(p)
}

make_type_1_error_med_h2_se_barplot <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=type_1_error)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=type_1_error_lb, ymax=type_1_error_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Type 1 Error") +
  		theme(legend.position="bottom")  +
  		geom_hline(yintercept=.05, linetype=2)

  	return(p)
}

 make_avg_h2_se_barplot <- function(df) {

 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000))
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
	 df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=med_h2)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Fraction h2 mediated\nby gene expression") +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=.1, linetype=2)
  	return(p)
}

make_avg_per_tissue_h2_se_barplot_for_single_class <- function(df, titler, sim_value, min_height, max_height) {
	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000))
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

#######################
# Command line args
#######################
global_simulation_name_string = args[1]
simulated_organized_results_dir = args[2]
visualize_simulated_results_dir = args[3]


if (FALSE) {
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

}

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


if (FALSE) {
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
}




