#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=75GB                         # Memory total in MiB (for all cores)



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
gene_trait_architecture="${20}"
eqtl_architecture="${21}"
eqtl_sample_size="${22}"
simulated_focus_input_dir="${23}"
simulated_focus_results_dir="${24}"
simulated_causal_twas_input_data_dir="${25}"
simulated_causal_twas_gene_models_dir="${26}"
simulated_causal_twas_results_dir="${27}"
processed_ctwas_genotype_data_dir="${28}"

source /home/bes710/.bash_profile
module load R/4.0.1
echo "Simulation"$simulation_number
date
echo $simulation_name_string
echo $eqtl_architecture
echo $eqtl_sample_size



# Make learned gene models output root seperate for each simulatino
simulated_learned_gene_models_dir=${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}"/"

global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"


##################################
# cTWAS
##################################
date
echo "Part 10: cTWAS"
sh run_single_causal_twas_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_causal_twas_input_data_dir $simulated_causal_twas_gene_models_dir $simulated_causal_twas_results_dir $eqtl_sample_size $processed_ctwas_genotype_data_dir

##################################
# FOCUS
##################################
date
echo "Part 11: FOCUS"
gene_type="component_gene"
sh run_single_focus_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir $gene_type








date
