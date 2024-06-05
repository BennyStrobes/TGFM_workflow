#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)




tmp_reference_genotype_stem="${1}"
simulated_gene_position_file="${2}"
processed_jlim_genotype_data_dir="${3}"
jlim_window_summary_file="${4}"



python3 get_reference_genotype_data_at_each_gene_for_jlim.py $tmp_reference_genotype_stem $simulated_gene_position_file $processed_jlim_genotype_data_dir $jlim_window_summary_file