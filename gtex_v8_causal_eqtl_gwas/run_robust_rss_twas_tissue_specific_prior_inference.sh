#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-40:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB   



trait_name="$1"
gtex_pseudotissue_file="$2"
pseudotissue_gtex_rss_multivariate_twas_dir="$3"
gene_version="$4"
gene_count_method="$5"
init_version="$6"

source ~/.bash_profile

python3 run_robust_rss_twas_tissue_specific_prior_inference.py $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_count_method $init_version

