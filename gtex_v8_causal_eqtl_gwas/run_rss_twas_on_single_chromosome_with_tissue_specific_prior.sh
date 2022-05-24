#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB   




chrom_num="$1"
trait_name="$2"
gtex_pseudotissue_file="$3"
component_data_file="$4"
ukbb_genome_wide_susie_organized_results_dir="$5"
pseudotissue_gtex_rss_multivariate_twas_dir="$6"


if false; then
source ~/.bash_profile
fi

echo $chrom_num
echo $trait_name

python3 run_rss_twas_on_single_chromosome_with_tissue_specific_prior.py $chrom_num $trait_name $gtex_pseudotissue_file $component_data_file $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir