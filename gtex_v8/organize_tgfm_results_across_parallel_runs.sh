#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=3GB                         # Memory total in MiB (for all cores)




tgfm_results_dir="$1"
gene_type="$2"
num_jobs="$3"
trait_names_file="$4"


source ~/.bash_profile

python3 organize_tgfm_results_across_parallel_runs.py $tgfm_results_dir $gene_type $num_jobs $trait_names_file