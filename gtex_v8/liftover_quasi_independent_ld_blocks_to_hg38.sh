#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)




liftover_directory="$1"
quasi_independent_dir="$2"
quasi_independent_ld_blocks_hg38_dir="$3"


source ~/.bash_profile

python3 liftover_quasi_independent_ld_blocks_to_hg38.py $liftover_directory $quasi_independent_dir $quasi_independent_ld_blocks_hg38_dir