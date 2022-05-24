#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB   



chrom_num="$1"
trait_name="$2"
ukbb_genome_wide_susie_organized_results_dir="$3"
gtex_pseudotissue_file="$4"
pseudotissue_gtex_susie_pmces_fusion_weights_dir="$5"
gtex_fusion_processed_intermediate_data="$6"
samp_size="$7"
pseudotissue_gtex_rss_multivariate_twas_data_dir="$8"
gene_version="$9"

echo $chrom_num"_"$trait_name"_"$gene_version

module load R/4.0.1

Rscript organize_data_for_univariate_rss_twas.R $chrom_num $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_rss_multivariate_twas_data_dir $gene_version


source ~/.bash_profile

python3 organize_data_for_multivariate_rss_twas.py $chrom_num $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_data_dir $samp_size $gene_version
