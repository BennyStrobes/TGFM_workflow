#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-9:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=60GB                         # Memory total in MiB (for all cores)



tgfm_input_summary_file="$1"
gtex_pseudotissue_file="$2"
output_file="$3"


source ~/.bash_profile
module load R/4.0.1

python3 compute_average_tissue_correlation_in_predicted_expression_across_genes.py $tgfm_input_summary_file $gtex_pseudotissue_file $output_file