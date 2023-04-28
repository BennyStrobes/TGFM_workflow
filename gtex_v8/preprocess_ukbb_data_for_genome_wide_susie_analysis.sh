#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=200GB                        # Memory total in MiB (for all cores)

source ~/.bash_profile

chrom_num="$1"
genome_wide_window_file="$2"
ukbb_sumstats_hg38_dir="$3"
gtex_genotype_dir="$4"
ref_1kg_genotype_dir="$5"
ukbb_preprocessed_for_genome_wide_susie_dir="$6"
ukbb_in_sample_ld_dir="$7"
ukbb_in_sample_genotype_dir="$8"

echo $chrom_num
if false; then
python3 preprocess_ukbb_data_for_genome_wide_susie_analysis.py $chrom_num $genome_wide_window_file $ukbb_sumstats_hg38_dir $gtex_genotype_dir $ref_1kg_genotype_dir $ukbb_preprocessed_for_genome_wide_susie_dir $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir
fi

python3 preprocess_ukbb_data_for_genome_wide_susie_analysis_pt_2_to_fix_ld.py $chrom_num $genome_wide_window_file $ukbb_sumstats_hg38_dir $gtex_genotype_dir $ref_1kg_genotype_dir $ukbb_preprocessed_for_genome_wide_susie_dir $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir
