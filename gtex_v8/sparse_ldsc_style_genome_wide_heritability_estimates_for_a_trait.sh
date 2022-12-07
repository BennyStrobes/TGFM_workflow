#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)





trait_name="$1"
ukkbb_window_summary_file="$2"
tissue_name_file="$3"
preprocessed_tgfm_data_dir="$4"
tgfm_heritability_results_dir="$5"
learn_intercept="$6"
output_stem="$7"
gene_type="$8"

date

if false; then
source ~/.bash_profile
module load R/4.0.1
fi



python3 sparse_ldsc_style_genome_wide_heritability_estimates_for_a_trait.py $trait_name $ukkbb_window_summary_file $tissue_name_file $preprocessed_tgfm_data_dir $tgfm_heritability_results_dir $learn_intercept $output_stem $gene_type

date