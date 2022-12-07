#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=50G                         # Memory total in MiB (for all cores)




baseline_ld_dir="$1"
output_root="$2"

source ~/.bash_profile


python3 preprocess_data_for_standard_sldsc.py $baseline_ld_dir $output_root