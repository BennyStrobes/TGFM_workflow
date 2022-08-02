#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB  


trait_file="$1"
gtex_pseudotissue_file="$2"
pseudotissue_gtex_rss_multivariate_twas_dir="$3"
rss_multivariate_twas_visualization_dir="$4"
trait_name="$5"


Rscript visualize_rss_twas.R $trait_file $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $rss_multivariate_twas_visualization_dir
