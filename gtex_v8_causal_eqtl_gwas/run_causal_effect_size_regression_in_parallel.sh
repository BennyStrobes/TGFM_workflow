#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)





trait_name="$1"
gtex_gene_file="$2"
gtex_susie_results_dir="$3"
ukbb_susie_results_dir="$4"
gtex_pseudotissue_file="$5"
gtex_onto_trait_causal_effect_size_regression_dir="$6"


source ~/.bash_profile

module load R/4.0.1


trait_name="$1"
gtex_gene_file="$2"
gtex_susie_results_dir="$3"
ukbb_susie_results_dir="$4"
gtex_pseudotissue_file="$5"
gtex_onto_trait_causal_effect_size_regression_dir="$6"
job_number="$7"
total_jobs="$8"


Rscript run_causal_effect_size_regression.R $trait_name $gtex_gene_file $gtex_susie_results_dir $ukbb_susie_results_dir $gtex_pseudotissue_file $gtex_onto_trait_causal_effect_size_regression_dir $job_number $total_jobs