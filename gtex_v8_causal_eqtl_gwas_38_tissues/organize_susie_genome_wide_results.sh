#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-60:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)






window_file="$1"
ukbb_genome_wide_susie_results_dir="$2"
ukbb_genome_wide_susie_organized_results_dir="$3"

source ~/.bash_profile
module load R/4.0.1
Rscript organize_susie_genome_wide_results.R $window_file $ukbb_genome_wide_susie_results_dir $ukbb_genome_wide_susie_organized_results_dir


source ~/.bash_profile
python3 filter_overlapping_components.py $window_file $ukbb_genome_wide_susie_organized_results_dir
