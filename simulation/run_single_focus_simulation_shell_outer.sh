#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=60GB                         # Memory total in MiB (for all cores)


simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
simulated_gene_expression_dir="$8"
simulated_learned_gene_models_dir="$9"
simulated_trait_dir="${10}"
simulated_gwas_dir="${11}"
simulated_tgfm_input_data_dir="${12}"
simulated_focus_input_dir="${13}"
simulated_focus_results_dir="${14}"



eqtl_sample_size="300"
sh run_single_focus_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir

eqtl_sample_size="500"
sh run_single_focus_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir

eqtl_sample_size="1000"
sh run_single_focus_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $eqtl_sample_size $simulated_focus_input_dir $simulated_focus_results_dir
