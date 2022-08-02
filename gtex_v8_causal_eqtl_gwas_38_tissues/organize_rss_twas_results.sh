#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





trait_name="$1"
gtex_pseudotissue_file="$2"
pseudotissue_gtex_rss_multivariate_twas_dir="$3"
gene_version="$4"
fusion_weights="$5"

if false; then
source ~/.bash_profile
fi
python3 organize_rss_twas_results.py $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $fusion_weights