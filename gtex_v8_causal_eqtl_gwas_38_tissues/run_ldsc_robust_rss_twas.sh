#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-50:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=60GB   



trait_name="$1"
gtex_pseudotissue_file="$2"
pseudotissue_gtex_rss_multivariate_twas_data_dir="$3"
ukbb_genome_wide_susie_organized_results_dir="$4"
pseudotissue_gtex_rss_multivariate_twas_dir="$5"
gene_version="$6"
chrom_num="$7"


source ~/.bash_profile
module load R/4.0.1


python3 run_ldsc_robust_rss_twas.py $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_data_dir $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $chrom_num
