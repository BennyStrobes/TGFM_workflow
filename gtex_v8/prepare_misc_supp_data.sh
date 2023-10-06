#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=3GB                         # Memory total in MiB (for all cores)



supp_data_dir="$1"
trait_names_file="$2"
gtex_pseudotissue_file="$3"
gtex_covariate_dir="$4"
tgfm_organized_results_dir="$5"


python3 prepare_misc_supp_data.py $supp_data_dir $trait_names_file $gtex_pseudotissue_file $gtex_covariate_dir $tgfm_organized_results_dir