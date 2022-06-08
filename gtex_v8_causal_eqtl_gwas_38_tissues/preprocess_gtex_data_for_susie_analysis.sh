#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-13:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=80G                         # Memory total in MiB (for all cores)


source ~/.bash_profile

chrom_num="$1"
gtex_tissue_file="$2"
xt_gene_list_file="$3"
eqtl_summary_stats_dir="$4"
ref_1kg_genotype_dir="$5"
ukbb_sumstats_hg38_dir="$6"
gtex_preprocessed_for_susie_dir="$7"

python3 preprocess_gtex_data_for_susie_analysis.py $chrom_num $gtex_tissue_file $xt_gene_list_file $eqtl_summary_stats_dir $ref_1kg_genotype_dir $ukbb_sumstats_hg38_dir $gtex_preprocessed_for_susie_dir