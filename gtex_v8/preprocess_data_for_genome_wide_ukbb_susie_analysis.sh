#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=160G                         # Memory total in MiB (for all cores)


source ~/.bash_profile

ukbb_sumstats_hg38_dir="$1"
gtex_genotype_dir="$2"
ref_1kg_genotype_dir="$3"
ukbb_preprocessed_for_genome_wide_susie_dir="$4"
ukbb_sumstats_hg19_dir="$5"
ukbb_in_sample_ld_dir="$6"
ukbb_in_sample_genotype_dir="$7"



##################################################
# Create genome-wide windows (should generate 2763) windows
# Actually less than 2763 because we are filtering out those in long LD regions
##################################################
genome_wide_window_file=$ukbb_preprocessed_for_genome_wide_susie_dir"genome_wide_susie_windows.txt"
if false; then
python3 generate_genome_wide_susie_windows.py $ukbb_sumstats_hg38_dir $genome_wide_window_file
fi

##################################################
# Generate trait file
##################################################
if false; then
python3 generate_trait_list_file_with_sample_size_and_heritabilities.py $ukbb_sumstats_hg38_dir $ukbb_sumstats_hg19_dir
fi

##################################################
# Generate processed input data for SuSiE for each window
##################################################
if false; then
for chrom_num in $(seq 1 22); do 
	sbatch preprocess_ukbb_data_for_genome_wide_susie_analysis.sh $chrom_num $genome_wide_window_file $ukbb_sumstats_hg38_dir $gtex_genotype_dir $ref_1kg_genotype_dir $ukbb_preprocessed_for_genome_wide_susie_dir $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir
done
fi


##################################################
# Concatenate window summary files across chromosomes
##################################################
if false; then
python3 merge_susie_input_window_file_across_chromosomes.py $ukbb_preprocessed_for_genome_wide_susie_dir
fi
