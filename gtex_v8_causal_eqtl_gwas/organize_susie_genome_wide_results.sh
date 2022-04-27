#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1






window_file="$1"
ukbb_genome_wide_susie_results_dir="$2"
ukbb_genome_wide_susie_organized_results_dir="$3"

if false; then
source ~/.bash_profile

module load R/4.0.1
Rscript organize_susie_genome_wide_results.R $window_file $ukbb_genome_wide_susie_results_dir $ukbb_genome_wide_susie_organized_results_dir
fi



source ~/.bash_profile
python3 check_for_overlapping_components_between_windows.py $window_file $ukbb_genome_wide_susie_organized_results_dir
