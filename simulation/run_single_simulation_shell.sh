#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=8GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
ldsc_real_data_results_dir="$8"
per_element_heritability="$9"
total_heritability="${10}"
fraction_expression_mediated_heritability="${11}"
ge_h2="${12}"
simulated_gene_expression_dir="${13}"
simulated_learned_gene_models_base_dir="${14}"
simulated_trait_dir="${15}"
simulated_gwas_dir="${16}"
simulated_tgfm_input_data_dir="${17}"
simulated_tgfm_results_dir="${18}"
simulated_coloc_results_dir="${19}"

source ~/.bash_profile
module load R/4.0.1
echo "Simulation"$simulation_number
date



mkdir ${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}
simulated_learned_gene_models_dir=${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}"/"

#######################################################
# Step 1: Simulate gene expression and fit gene models
#######################################################
echo "Simulation Step 1"
python3 simulate_gene_expression_and_fit_gene_model.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $ge_h2



#######################################################
# Step 2: Simulate trait values
#######################################################
echo "Simulation Step 2"
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $simulated_trait_dir $n_gwas_individuals



source ~/.bash_profile
#######################################################
# Step 8: Run GWAS on simulated trait on only snps in TGFM windows.
#######################################################
echo "Simulation Step 8"
global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"
python3 run_gwas_on_simulated_trait_at_snps_in_tgfm_windows.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $global_window_file $simulated_gwas_dir


# Merge gwas data across windows
source ~/.bash_profile
merged_gwas_summary_stat_file=${simulated_gwas_dir}${simulation_name_string}"_merged_gwas_summary_stats.txt"
python3 generate_merged_gwas_data.py $global_window_file $simulation_number $chrom_num $simulation_name_string ${simulated_gwas_dir} $processed_genotype_data_dir $n_gwas_individuals $merged_gwas_summary_stat_file


#######################################################
# Step 9: Run coloc
#######################################################
echo "Simulation Step 9"
source ~/.bash_profile
module load R/4.0.1
eqtl_sample_size_arr=( "300" "500" "1000")
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 run_coloc_shell.py $merged_gwas_summary_stat_file $simulation_number $chrom_num $simulation_name_string $eqtl_sample_size $n_gwas_individuals $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_coloc_results_dir
done



source ~/.bash_profile

#######################################################
# Step 10: Preprocess data for TGFM
#######################################################
echo "Simulation Step 10"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_type="susie"



eqtl_sample_size_arr=( "300" "500" "1000")
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	echo $eqtl_sample_size
	python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $global_window_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done


# Can probabably delete eqtl susie outputs
# Can delete per window summary stats (will be saved elsewhere)


date



