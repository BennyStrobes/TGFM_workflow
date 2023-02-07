#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




chrom_num="$1"
cis_window="$2"
n_gwas_individuals="$3"
global_simulation_name_string="$4"
total_heritability="$5"
fraction_expression_mediated_heritability="$6"
simulated_sldsc_results_dir="$7"
simulated_organized_results_dir="$8"
simulated_tgfm_results_dir="$9"
simulated_trait_dir="${10}"


if false; then
source ~/.bash_profile
fi


python3 organize_simulation_results_across_parallel_simulations.py $chrom_num $cis_window $n_gwas_individuals $global_simulation_name_string $total_heritability $fraction_expression_mediated_heritability $simulated_sldsc_results_dir $simulated_organized_results_dir $simulated_tgfm_results_dir $simulated_trait_dir