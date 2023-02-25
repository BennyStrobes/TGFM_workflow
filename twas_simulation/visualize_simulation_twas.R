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






make_log_bayes_factor_scatterplot <- function(input_file, model_name, eqtl_sample_size) {
	df <- read.table(input_file, header=TRUE)
	df$pmces_log_p_value = -log(df$pmces_p_value)
	df$distribution_log_p_value = -log(df$distribution_p_value)

	p <- ggplot(df, aes(x=distribution_log_p_value, y=pmces_log_p_value, color=variance_ratio)) + geom_point() +
	  figure_theme() +
	  labs(x="Gene distribution -log10(p)", y="Gene PMCES -log10(p)", fill="Expression variance\nratio", title=paste0(model_name, " / ", eqtl_sample_size))
	return(p)
}



make_power_fdr_line_plot <- function(input_file, eqtl_sample_size) {
	df <- read.table(input_file, header=TRUE)
	df = df[as.character(df$eqtl_sample_size)==eqtl_sample_size,]
	p <- ggplot(data=df, aes(x=fdr, y=power, group=twas_model)) +
  		geom_line(aes(color=twas_model)) +
  		figure_theme() +
  		labs(title="eQTL sample size")

  	return(p)
}

make_fdr_barplot <- function(fdr_input_file, threshold) {
	df <- read.table(fdr_input_file, header=TRUE)
	df$eqtl_sample_size = factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000))
	df = df[df$eqtl_sample_size %in% c(100,200,300,500,1000),]
	df$twas_model = factor(df$twas_model, levels=c("marginal_pmces", "marginal_distr", "susie_pmces", "susie_distr", "fusion_lasso_pmces"))

	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=fdr, fill=twas_model)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=fdr_lb, ymax=fdr_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="FDR", fill="")  +
  		geom_hline(yintercept=(threshold), linetype=2) +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)


}


make_power_barplot <- function(fdr_input_file, threshold) {
	df <- read.table(fdr_input_file, header=TRUE)
	df$eqtl_sample_size = factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000))
	df = df[df$eqtl_sample_size %in% c(100,200,300,500,1000),]
	df$twas_model = factor(df$twas_model, levels=c("marginal_pmces", "marginal_distr", "susie_pmces", "susie_distr", "fusion_lasso_pmces"))

	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=power, fill=twas_model)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="Power", fill="")  +
  		theme(plot.title = element_text(hjust = 0.5,size=12)) +
  		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)


}


simulation_name_string = args[1]
simulated_organized_results_dir = args[2]
visualize_twas_results_dir = args[3]


#####################
# FDR plot se bar plot
#####################
p_value_thresh="0.05"
fdr_input_file <- paste0(simulated_organized_results_dir, simulation_name_string, "_statistical_fdr_p_", p_value_thresh, ".txt")
fdr_barplot <- make_fdr_barplot(fdr_input_file, 0.05)
output_file <- paste0(visualize_twas_results_dir, "fdr_barplot_", p_value_thresh,".pdf")
ggsave(fdr_barplot, file=output_file, width=7.2, height=3.7, units="in")

#####################
# power plot se bar plot
#####################
p_value_thresh="0.05"
power_input_file <- paste0(simulated_organized_results_dir, simulation_name_string, "_statistical_power_p_", p_value_thresh, ".txt")
power_barplot <- make_power_barplot(power_input_file, 0.05)
output_file <- paste0(visualize_twas_results_dir, "power_barplot_", p_value_thresh,".pdf")
ggsave(power_barplot, file=output_file, width=7.2, height=3.7, units="in")




#####################
# Power-FDR plot at various eQTL sample size
#####################
eqtl_sample_sizes <- c("100", "200", "300", "500", "1000")
for (eqtl_sample_iter in 1:length(eqtl_sample_sizes)) {
	eqtl_sample_size = eqtl_sample_sizes[eqtl_sample_iter]

	input_file <- paste0(simulated_organized_results_dir, simulation_name_string, "_power_false_discovery_curve_input.txt")

	power_fdr_scatter <- make_power_fdr_line_plot(input_file, eqtl_sample_size)
	output_file <- paste0(visualize_twas_results_dir, "power_fdr_curve_", eqtl_sample_size, ".pdf")
	ggsave(power_fdr_scatter, file=output_file, width=7.2, height=5.7, units="in")

}



#####################
# Make log bayes factor plot
#####################
eqtl_sample_sizes <- c("100", "200", "300", "500", "1000")
model_names <- c("susie", "marginal")
for (eqtl_sample_iter in 1:length(eqtl_sample_sizes)) {
	eqtl_sample_size = eqtl_sample_sizes[eqtl_sample_iter]
	for (model_iter in 1:length(model_names)) {
		model_name <- model_names[model_iter]


		input_file <- paste0(simulated_organized_results_dir, simulation_name_string, "_organized_bf_comparison_", model_name, "_", eqtl_sample_size, ".txt")
		# make log bayes factor plot
		bf_scatter <- make_log_bayes_factor_scatterplot(input_file, model_name, eqtl_sample_size)
		output_file <- paste0(visualize_twas_results_dir, "log_bf_scatter_comparison_", model_name, "_", eqtl_sample_size, ".pdf")
		ggsave(bf_scatter, file=output_file, width=7.2, height=5.7, units="in")

	}
}