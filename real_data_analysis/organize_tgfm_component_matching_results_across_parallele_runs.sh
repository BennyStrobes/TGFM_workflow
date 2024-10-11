#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-02:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)


tgfm_results_dir="$1"
gene_type="$2"
num_jobs="$3"
trait_names_file="$4"
tgfm_organized_results_dir="$5"



python3 organize_tgfm_component_matching_results_across_parallele_runs.py $tgfm_results_dir $gene_type $num_jobs $trait_names_file $tgfm_organized_results_dir