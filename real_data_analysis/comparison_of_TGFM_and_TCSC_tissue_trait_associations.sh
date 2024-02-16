#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4GB                         # Memory total in MiB (for all cores)





tgfm_results_dir="$1"
trait_names_file="$2"
tcsc_results_file="$3"
tcsc_comparison_dir="$4"


python3 comparison_of_TGFM_and_TCSC_tissue_trait_associations.py $tgfm_results_dir $trait_names_file $tcsc_results_file $tcsc_comparison_dir