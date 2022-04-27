#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-13:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=80GB                         # Memory total in MiB (for all cores)



gtex_gene_file="$1"
ukbb_sumstats_hg38_dir="$2"
ukbb_preprocessed_for_susie_dir="$3"


source ~/.bash_profile

python3 preprocess_ukbb_data_for_susie_analysis.py $gtex_gene_file $ukbb_sumstats_hg38_dir $ukbb_preprocessed_for_susie_dir