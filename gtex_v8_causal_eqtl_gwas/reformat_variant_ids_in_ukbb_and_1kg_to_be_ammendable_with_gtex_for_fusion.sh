#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)



gtex_genotype_dir="$1"
ukbb_sumstats_hg38_dir="$2"
ref_1kg_genotype_dir="$3"
gtex_fusion_processed_intermediate_data="$4"


source ~/.bash_profile
if false; then
python3 reformat_ukbb_sumstats_to_be_ammendable_with_gtex_for_fusion.py $gtex_genotype_dir $ukbb_sumstats_hg38_dir $gtex_fusion_processed_intermediate_data
fi
if false; then
python3 reformat_1kg_genotypes_to_be_ammendable_with_gtex_for_fusion.py $gtex_genotype_dir $ref_1kg_genotype_dir $gtex_fusion_processed_intermediate_data
fi