#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-9:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=14GB                         # Memory total in MiB (for all cores)



simulated_twas_dir="$1"
simulated_organized_results_dir="$2"
simulation_name_string="$3"
simulated_expr_snp_corr_dir="$4"




python3 organize_simulated_twas_results.py $simulated_twas_dir $simulated_organized_results_dir $simulation_name_string $simulated_expr_snp_corr_dir
