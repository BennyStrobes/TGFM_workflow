#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)


num_jobs="$1"
gene_type="$2"
preprocessed_tgfm_data_dir="$3"


source ~/.bash_profile


python3 organize_processed_tgfm_input_data.py $num_jobs $gene_type $preprocessed_tgfm_data_dir
