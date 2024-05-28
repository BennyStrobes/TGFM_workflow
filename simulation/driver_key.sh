#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)


#############################
# Simulation of:
####### 1. TGFM-SLDSC
####### 2. TGFM fine-mapping
##############################




############################
# Input data
############################

# Directory containing ldsc scripts
ldsc_code_dir="/n/groups/price/ben/tools/ldsc/"

# Mod-ldscore regression code
curr_dir=`pwd`
mod_ldsc_code_dir=$curr_dir"/modified_sldsc/"

# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"


# Ldsc weights for hg19
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/"

# LDSC baseline annotations (hg19)
ldsc_baseline_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline_v1.2/"

# LDSC baselineLD annotations (hg19)
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# LDSC 1KG genotype files (hg19)
kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"

# Directory containing LDSC Summary statistics (This will be used to create a realistic gene model)
ldsc_summary_stats_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"



############################
# Output directories
############################
# Output roots (one temporary on scratch and one permanent)
temp_output_root="/n/scratch/users/b/bes710/causal_eqtl_gwas/simulation/"
alkes_temp_output_root="/n/scratch/users/a/ap92/ben/tgfm/simulation/"
perm_output_root="/n/groups/price/ben/causal_eqtl_gwas/simulation/"

# Directory containing processed genotype data
ldsc_real_data_results_dir=$temp_output_root"ldsc_real_data_results/"

# Directory containing processed genotype data
processed_genotype_data_dir=$temp_output_root"processed_genotype/"

# Directory containing processed genotype data
processed_ctwas_genotype_data_dir=$temp_output_root"processed_ctwas_genotype/"

# Directory containing simulated gene positions
simulated_gene_position_dir=$temp_output_root"simulated_gene_positions/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_gene_expression_dir=$perm_output_root"simulated_gene_expression/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_learned_gene_models_dir=$perm_output_root"simulated_learned_gene_models/"

# Directory containing simulated trait values
simulated_trait_dir=$temp_output_root"simulated_trait/"

# Directory contaiing simulated gwas results
simulated_gwas_dir=$temp_output_root"simulated_gwas/"

# Directory containing simulated coloc results
simulated_coloc_results_dir=$temp_output_root"simulated_coloc_results/"

# Directory containing organized simulation results
simulated_organized_results_dir=$perm_output_root"simulated_organized_results/"

# Directory containing simulated tgfm input data
simulated_tgfm_input_data_dir=$alkes_temp_output_root"simulated_tgfm_input/"

# Directory containing simulated tgfm results
simulated_tgfm_results_dir=$alkes_temp_output_root"simulated_tgfm_results/"

# Directory containing simulated tgfm-two-step results
simulated_two_step_tgfm_results_dir=$alkes_temp_output_root"simulated_tgfm_two_step_results/"

# Directory containing simulated best tagging gene-tissue pairs
simulated_best_tagging_gt_dir=$alkes_temp_output_root"simulated_best_tagging_gene_tissue_pairs/"

# Directory containing simulated focus input data
simulated_focus_input_dir=$temp_output_root"simulated_focus_input/"

# Directory containing simulated focus results dir
simulated_focus_results_dir=$temp_output_root"simulated_focus_results/"

# Directory containing simulated causal-twas input data
simulated_causal_twas_input_data_dir=$temp_output_root"simulated_causal_twas_input/"

# Directory containing simulated causal-twas gene models data
simulated_causal_twas_gene_models_dir=$perm_output_root"simulated_causal_twas_gene_models/"

# Directory containing simulated causal-twas results
simulated_causal_twas_results_dir=$temp_output_root"simulated_causal_twas_results/"

# Directory containing visualizations of simulated results
visualize_simulated_results_dir=$perm_output_root"visualize_simulated_results/"
visualize_simulated_results_debug_dir=$perm_output_root"visualize_simulated_results_debug/"



############################
# Simulation parameters
############################
# Chromosome to simulate on 
chrom_num="1"





############################
# Run S-LDSC on real data.
## This is done in order to have a realistic heritability model (taus)
############################
if false; then
sh run_ldsc_on_real_trait_data_to_get_realistic_heritability_model.sh $ldsc_code_dir $ldsc_baseline_hg19_annotation_dir $ldsc_weights_dir $ldsc_summary_stats_dir $kg_genotype_dir $ldsc_real_data_results_dir
fi

############################
# Prepare genotype data for analysis:
## 1. Filter number of individuals in original data
## 2. Filter sites to be those in LDSC annotation file
## 3. Convert to plink bed files
############################
# NOTE: THERE IS CURRENTLY A HACK IN HERE TO REMOVE 3 variants (out of 500000) on chrom 1 that have no variance across the 100-sample eqtl data set.
############################
n_gwas_individuals="100000"
if false; then
sbatch prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir $processed_ctwas_genotype_data_dir
fi

if false; then
n_gwas_individuals="200000"
sbatch prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir $processed_ctwas_genotype_data_dir

n_gwas_individuals="50000"
sbatch prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir $processed_ctwas_genotype_data_dir
fi

############################
# Prepare gene file for simulation:
# Genes are defined by actual tss
# Limit to protein coding genes
# In simulation, I will assume gene is expressed in each tissue
############################
simulated_gene_position_file=${simulated_gene_position_dir}"gene_positions_chr"${chrom_num}".txt"
if false; then
sh prepare_simulated_gene_position_list.sh $chrom_num $gencode_gene_annotation_file $simulated_gene_position_file
fi




############################
# Main simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# cis window arround genes to define eQTLs
cis_window="100000"

# Gene expression heritability
ge_h2="075"

# Gene-trait architecture
gene_trait_architecture="1_caus_t"
gene_trait_architecture="2_caus_t"

# eQTL architecture
eqtl_architecture="random_Neqtl"
eqtl_architecture="selection_1"
eqtl_architecture="selection_125"
eqtl_architecture="default"


############################
# Run main single simulation of simulating the trait data
############################
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi


############################
# Fit gene models in single fine-mapping simulation
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "realistic" "100" "300" "500" "1000")
 # Done seperately due to drastically diff time limits
eqtl_sample_size_arr=( "500" "1000") # 60 h
eqtl_sample_size_arr=( "realistic" "100" "300" ) # 20h
run_lasso="run_lasso"
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_learn_gene_models_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $run_lasso
done
done
fi




############################
# Run main single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "realistic" "100" "300" "500" "1000")
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi



############################
# Run method-comparison single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
simulation_number="1"
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "realistic" "100" "300" "500" "1000")
if false; then
for simulation_number in $(seq 6 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_method_comparison_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi


############################
# Run comparison to two-step fine-mapping approach
############################
simulation_number="1"
eqtl_sample_size_arr=( "realistic" "100" "300" "500" "1000")
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_two_step_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_two_step_tgfm_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size
done 
done
fi




############################
# Extra run of creating gene models for varying number of samples simulation
############################
eqtl_sample_size="500"
run_lasso="no_run_lasso"  # reduces running time to like 15 h
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_learn_gene_models_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $run_lasso
done
fi

############################
# Run single TGFM fine-mapping simulation for varying number of samples simulation
############################
eqtl_sample_size="500"
if false; then
for simulation_number in $(seq 1 100); do 
	n_samples="50"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $n_samples

	n_samples="200"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $n_samples
done
fi

############################
# Run single TGFM fine-mapping simulation for varying the prior
############################
eqtl_sample_size="500"
n_samples="100"
if false; then
for simulation_number in $(seq 1 50); do 
	nm_var_prior_multiplier="0.5"
	gene_tissue_prior_multiplier="1.0"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_fine_mapping_simulation_altering_prior.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $n_samples $nm_var_prior_multiplier $gene_tissue_prior_multiplier

	nm_var_prior_multiplier="2.0"
	gene_tissue_prior_multiplier="1.0"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_fine_mapping_simulation_altering_prior.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $n_samples $nm_var_prior_multiplier $gene_tissue_prior_multiplier

	nm_var_prior_multiplier="1.0"
	gene_tissue_prior_multiplier="0.5"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_fine_mapping_simulation_altering_prior.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $n_samples $nm_var_prior_multiplier $gene_tissue_prior_multiplier

	nm_var_prior_multiplier="1.0"
	gene_tissue_prior_multiplier="2.0"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_fine_mapping_simulation_altering_prior.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $n_samples $nm_var_prior_multiplier $gene_tissue_prior_multiplier
done
fi




############################
# Selection simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# cis window arround genes to define eQTLs
cis_window="100000"

# Gene expression heritability
ge_h2="075"

# Gene-trait architecture
gene_trait_architecture="1_caus_t"
gene_trait_architecture="2_caus_t"

# eQTL architecture
eqtl_architecture="random_Neqtl"
eqtl_architecture="selection_1"


############################
# Run simulating the trait data
############################
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi



############################
# Fit gene models in single fine-mapping simulation
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "300" "500") # 45 h
run_lasso="run_lasso"
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_learn_gene_models_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $run_lasso
done
done
fi


############################
# Run main single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "300" "500")
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi

if false; then
sed 1d "/n/groups/price/ben/code_temp7/missing_runs.txt" | while read simulation_number eqtl_sample_size; do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
fi


############################
# Run method-comparison single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
simulation_number="1"
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "300" "500") 
if false; then
for simulation_number in $(seq 1 100); do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_method_comparison_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi




############################
# Pleiotropy  simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# cis window arround genes to define eQTLs
cis_window="100000"

# Gene expression heritability
ge_h2="075"

# Gene-trait architecture
gene_trait_architecture="1_caus_t"
gene_trait_architecture="2_caus_t"

# eQTL architecture
eqtl_architecture="random_Neqtl"
eqtl_architecture="selection_1"
eqtl_architecture="pleiotropy"


############################
# Run simulating the trait data
############################
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi


############################
# Fit gene models in single fine-mapping simulation
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "300" "500") # 45 h
run_lasso="no_run_lasso"
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_learn_gene_models_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $run_lasso
done
done
fi


############################
# Run main single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "300" "500")
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi


############################
# Random_Neqtl simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# cis window arround genes to define eQTLs
cis_window="100000"

# Gene expression heritability
ge_h2="075"

# Gene-trait architecture
gene_trait_architecture="1_caus_t"
gene_trait_architecture="2_caus_t"

# eQTL architecture
eqtl_architecture="selection_1"
eqtl_architecture="pleiotropy"
eqtl_architecture="random_Neqtl"


############################
# Run simulating the trait data
############################
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi


############################
# Fit gene models in single fine-mapping simulation
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "300" "500") # 45 h
run_lasso="no_run_lasso"
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_learn_gene_models_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $run_lasso
done
done
fi



############################
# Run main single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "300" "500")
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi





############################
# Main simulation with 1 causal tissue simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# cis window arround genes to define eQTLs
cis_window="100000"

# Gene expression heritability
ge_h2="075"

# Gene-trait architecture
gene_trait_architecture="2_caus_t"
gene_trait_architecture="1_caus_t"

# eQTL architecture
eqtl_architecture="random_Neqtl"
eqtl_architecture="selection_1"
eqtl_architecture="selection_125"
eqtl_architecture="default"


############################
# Run main single simulation of simulating the trait data
############################
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi


############################
# Fit gene models in single fine-mapping simulation
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "realistic" "100" "300" "500" "1000")
 # Done seperately due to drastically diff time limits
eqtl_sample_size_arr=( "500" "1000") # 60 h
eqtl_sample_size_arr=( "realistic" "300" ) # 20h

run_lasso="run_lasso"
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_learn_gene_models_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir $run_lasso
done
done
fi

############################
# Run main single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "realistic" "300" "500" "1000")
if false; then
for simulation_number in $(seq 1 100); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi

############################
# Run method-comparison single fine-mapping simulation of the simulated trait data at a single eQTL sample
############################
# Loop through eqtl sample sizes
eqtl_sample_size_arr=( "realistic" "300" "500" "1000")
if false; then
for simulation_number in $(seq 1 100); do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_method_comparison_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $processed_ctwas_genotype_data_dir
done
done
fi

############################
# Run comparison to two-step fine-mapping approach
############################
simulation_number="1"
eqtl_sample_size_arr=( "realistic" "300" "500" "1000")
if false; then
for simulation_number in $(seq 1 1); do 
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_two_step_fine_mapping_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_two_step_tgfm_results_dir $gene_trait_architecture $eqtl_architecture $eqtl_sample_size
done 
done
fi













############################
# Organize simulation results across parallel simulations
############################
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
if false; then
sh organize_simulation_results_across_parallel_simulations.sh $chrom_num $cis_window $n_gwas_individuals $global_simulation_name_string $total_heritability $fraction_expression_mediated_heritability $simulated_organized_results_dir $simulated_tgfm_results_dir $simulated_trait_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $simulated_gene_position_file $processed_genotype_data_dir $simulated_focus_results_dir $simulated_coloc_results_dir $simulated_causal_twas_results_dir $simulated_two_step_tgfm_results_dir
fi







############################
# Visualize single simulation
############################
if false; then
source /home/bes710/.bash_profile
module load R/3.5.1
fi
if false; then
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
Rscript visualize_single_simulation.R $global_simulation_name_string $simulated_organized_results_dir $visualize_simulated_results_dir
fi




























###########
# OLD
###########

############################
# Run main single simulation of processing
############################
# Iteration of simulation (also works as seed)
eqtl_architecture="default"
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi




if false; then
for simulation_number in $(seq 101 200); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi

if false; then
for simulation_number in $(seq 201 210); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi

if false; then
for simulation_number in $(seq 211 300); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi

if false; then
	simulation_number="81"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
	
	simulation_number="60"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
fi


eqtl_architecture="random_Neqtl"
if false; then
for simulation_number in $(seq 1 50); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi



eqtl_architecture="selection_1"
if false; then
for simulation_number in $(seq 1 50); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi


eqtl_architecture="selection_125"
if false; then
for simulation_number in $(seq 1 50); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi



############################
# Run single simulation of TGFM
############################
tgfm_tissues="all_t"

eqtl_architecture="default"

if false; then
for simulation_number in $(seq 1 100); do 
	eqtl_sample_size="100"
	gene_type="component_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi

if false; then
for simulation_number in $(seq 2 100); do 
	eqtl_sample_size="100"
	gene_type="all_non_zero_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi

if false; then
eqtl_sample_size="100"
for simulation_number in $(seq 1 59); do 
	gene_type="all_non_zero_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type

	gene_type="component_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done

for simulation_number in $(seq 61 80); do 
	gene_type="all_non_zero_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type

	gene_type="component_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done

for simulation_number in $(seq 82 100); do 
	gene_type="all_non_zero_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type

	gene_type="component_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done


eqtl_sample_size="1000"
for simulation_number in $(seq 1 59); do 
	gene_type="all_non_zero_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type

	gene_type="component_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done

for simulation_number in $(seq 61 80); do 
	gene_type="all_non_zero_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type

	gene_type="component_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done

for simulation_number in $(seq 82 100); do 
	gene_type="all_non_zero_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type

	gene_type="component_gene"
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi





eqtl_sample_size="100"
gene_type="max_min_ratio_5"
if false; then
for simulation_number in $(seq 101 200); do
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi




gene_type="max_min_ratio_50"
if false; then
for simulation_number in $(seq 101 200); do
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi


eqtl_sample_size="100"
gene_type="cafeh"
if false; then
for simulation_number in $(seq 201 300); do
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_w_cafeh_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi




if false; then
for simulation_number in $(seq 2 50); do 
	# Simulation string used for output file

	eqtl_sample_size="300"
	eqtl_architecture="random_Neqtl"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir

	eqtl_architecture="selection_1"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir

	eqtl_architecture="selection_125"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir


	eqtl_sample_size="500"
	eqtl_architecture="random_Neqtl"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir

	eqtl_architecture="selection_1"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir

	eqtl_architecture="selection_125"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
done
fi

if false; then
for simulation_number in $(seq 2 17); do 
	eqtl_sample_size="500"
	eqtl_architecture="selection_1"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell_all_genes.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
done

for simulation_number in $(seq 19 50); do 
	eqtl_sample_size="500"
	eqtl_architecture="selection_1"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell_all_genes.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
done
fi

if false; then
for simulation_number in $(seq 1 1); do 
	eqtl_sample_size="500"
	eqtl_architecture="selection_1"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell_all_genes.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
done
fi

if false; then
for simulation_number in $(seq 18 18); do 
	eqtl_sample_size="500"
	eqtl_architecture="selection_1"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell_all_genes.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
done
fi

############################
# Run single simulation of TGFM
############################
# cis window arround genes to define eQTLs
cis_window="100000"

if false; then
tgfm_tissues="all_t"
tgfm_tissues="no_t0"
simulation_number="2"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}
eqtl_sample_size="300"
sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
fi


if false; then
for simulation_number in $(seq 31 100); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}

	eqtl_sample_size="300"
	tgfm_tissues="all_t"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
	tgfm_tissues="no_t0"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
		
	eqtl_sample_size="500"
	tgfm_tissues="all_t"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
	tgfm_tissues="no_t0"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir

	eqtl_sample_size="1000"
	tgfm_tissues="all_t"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
	tgfm_tissues="no_t0"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
done
fi



if false; then
simulation_number="25"
eqtl_sample_size="1000"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}
tgfm_tissues="no_t0"
sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size $tgfm_tissues $simulated_best_tagging_gt_dir
fi





# SIMULATION PREPROCESSING with realistic eqtl sample sizes
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"
# cis window arround genes to define eQTLs
cis_window="100000"
# Per genetic-element heritabilities
per_element_heritability="0.0005"
# Total heritability
total_heritability="0.3"
# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"
# cis window arround genes to define eQTLs
cis_window="100000"
# Gene expression heritability
ge_h2="075"
# Gene-trait architecture
gene_trait_architecture="1_caus_t"
# eQTL architecture
eqtl_architecture="default"

if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_realistic_qtl_ss_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_simulation_shell_realistic_qtl_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi
if false; then
for simulation_number in $(seq 101 200); do 
simulation_name_string="simulation_realistic_qtl_ss_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
sbatch run_single_simulation_shell_realistic_qtl_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir $gene_trait_architecture $eqtl_architecture
done
fi

############################
# Run single simulation of TGFM
############################
tgfm_tissues="all_t"
eqtl_architecture="default"
if false; then
for simulation_number in $(seq 1 100); do 
	gene_type="component_gene"
	simulation_name_string="simulation_realistic_qtl_ss_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell_realistic_qtl_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi

if false; then
for simulation_number in $(seq 101 200); do 
	gene_type="component_gene"
	simulation_name_string="simulation_realistic_qtl_ss_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_tgfm_simulation_shell_realistic_qtl_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $tgfm_tissues $simulated_best_tagging_gt_dir $gene_type
done
fi




















# Alt parameters
if false; then
n_gwas_individuals="50000"
for simulation_number in $(seq 2 50); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir
done

n_gwas_individuals="200000"
for simulation_number in $(seq 2 50); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir
done

ge_h2="05"
n_gwas_individuals="100000"
for simulation_number in $(seq 2 50); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir
done

ge_h2="1"
n_gwas_individuals="100000"
for simulation_number in $(seq 2 50); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir
done
fi



############################
# Run single simulation of TGFM
############################
# cis window arround genes to define eQTLs
cis_window="100000"

if false; then
for simulation_number in $(seq 1 20); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}

	eqtl_sample_size="300"
	sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
	
	eqtl_sample_size="500"
	sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size

	eqtl_sample_size="1000"
	sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
done
fi


if false; then
for simulation_number in $(seq 21 100); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}

	eqtl_sample_size="300"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
	
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size

	eqtl_sample_size="1000"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
done
fi




if false; then
for simulation_number in $(seq 1 20); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
done
fi



############################
# Run single simulation of TGFM with alt. parameters
############################
# cis window arround genes to define eQTLs
cis_window="100000"

if false; then
for simulation_number in $(seq 1 50); do 
	n_gwas_individuals="50000"
	ge_h2="075"
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size

	n_gwas_individuals="200000"
	ge_h2="075"
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size

	ge_h2="05"
	n_gwas_individuals="100000"	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size

	ge_h2="1"
	n_gwas_individuals="100000"	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
done
fi

############################
# Run FOCUS simulation
############################
if false; then
for simulation_number in $(seq 1 20); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_focus_simulation_shell_outer.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_focus_input_dir $simulated_focus_results_dir
done
fi

if false; then
for simulation_number in $(seq 21 100); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_focus_simulation_shell_outer.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_focus_input_dir $simulated_focus_results_dir
done
fi

############################
# Run causal-TWAS simulation
############################
gene_trait_architecture="1_caus_t"
eqtl_architecture="default"
eqtl_sample_size="realistic_qtl_ss"
simulation_number="101"
simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
if false; then
sbatch run_single_causal_twas_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $eqtl_sample_size $processed_ctwas_genotype_data_dir
fi

if false; then
for simulation_number in $(seq 103 200); do 
	simulation_name_string="simulation_new_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_gt_arch_"${gene_trait_architecture}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_causal_twas_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $eqtl_sample_size $processed_ctwas_genotype_data_dir
done
fi



# Organize simulation results across parallel simulations
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
if false; then
sh organize_simulation_results_across_parallel_simulations.sh $chrom_num $cis_window $n_gwas_individuals $global_simulation_name_string $total_heritability $fraction_expression_mediated_heritability $simulated_organized_results_dir $simulated_tgfm_results_dir $simulated_trait_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $simulated_gene_position_file $processed_genotype_data_dir $simulated_focus_results_dir $simulated_coloc_results_dir $simulated_causal_twas_results_dir
fi




############################
# Visualize single simulation
############################
if false; then
source /home/bes710/.bash_profile
module load R/3.5.1
fi
if false; then
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
Rscript visualize_single_simulation.R $global_simulation_name_string $simulated_organized_results_dir $visualize_simulated_results_dir
fi

if false; then
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
Rscript visualize_single_simulation_debug.R $global_simulation_name_string $simulated_organized_results_dir $visualize_simulated_results_debug_dir
fi
