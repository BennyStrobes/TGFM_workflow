#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4GB                         # Memory total in MiB (for all cores)



trait_names_file="$1"
orig_tgfm_results_dir="$2"
new_tgfm_results_dir="$3"
original_tissue_names="$4"
new_tissue_names="$5"
removed_tissue_name="$6"
replication_output_root="$7"
suffix1="$8"
suffix2="$9"


python3 tgfm_tissue_replication_analysis.py $trait_names_file $orig_tgfm_results_dir $new_tgfm_results_dir $original_tissue_names $new_tissue_names $removed_tissue_name $replication_output_root $suffix1 $suffix2