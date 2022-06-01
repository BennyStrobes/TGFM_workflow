#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-50:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB   




chrom_num="$1"
trait_name="$2"
gtex_pseudotissue_file="$3"
component_data_file="$4"
ukbb_genome_wide_susie_organized_results_dir="$5"
pseudotissue_gtex_rss_multivariate_twas_dir="$6"
gene_version="$7"


source ~/.bash_profile
module load R/4.0.1


echo $chrom_num
echo $trait_name
echo $gene_version

python3 run_robust_rss_twas_on_single_chromosome.py $chrom_num $trait_name $gtex_pseudotissue_file $component_data_file $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version