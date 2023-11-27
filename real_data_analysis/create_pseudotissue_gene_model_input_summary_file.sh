#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4GB                         # Memory total in MiB (for all cores)


pseudotissue_name="$1"
composit_tissue_string="$2"
gtex_processed_expression_dir="$3"
gtex_pseudotissue_gene_model_input_dir="$4"
num_jobs="$5"

source ~/.bash_profile



python3 create_pseudotissue_gene_model_input_summary_file.py $pseudotissue_name $composit_tissue_string $gtex_processed_expression_dir $gtex_pseudotissue_gene_model_input_dir $num_jobs