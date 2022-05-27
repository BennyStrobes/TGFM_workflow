#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)





gtex_susie_input_data_file="$1"
trait_name="$2"
sumstat_file="$3"
trait_sample_size="$4"
coloc_data_dir="$5"

source ~/.bash_profile



python3 preprocess_coloc_data_for_single_trait.py $gtex_susie_input_data_file $trait_name $sumstat_file $trait_sample_size $coloc_data_dir