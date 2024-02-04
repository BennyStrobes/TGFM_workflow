#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-18:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=200GB                        # Memory total in MiB (for all cores)

source ~/.bash_profile

chrom_num="$1"
genome_wide_window_file="$2"
ukbb_sumstats_hg38_dir="$3"
ukbb_preprocessed_for_genome_wide_susie_dir="$4"
ukbb_in_sample_ld_dir="$5"
ukbb_in_sample_genotype_dir="$6"

echo $chrom_num
python3 preprocess_ukbb_data_for_genome_wide_susie_analysis.py $chrom_num $genome_wide_window_file $ukbb_sumstats_hg38_dir $ukbb_preprocessed_for_genome_wide_susie_dir $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir

