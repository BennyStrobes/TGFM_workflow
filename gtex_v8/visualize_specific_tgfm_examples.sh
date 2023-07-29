#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=15GB                         # Memory total in MiB (for all cores)

specific_examples_input_file="$1"
tgfm_input_summary_file="$2"
tgfm_results_dir="$3"
tgfm_organized_results_dir="$4"
gtex_susie_gene_models_dir="$5"
gene_annotation_file="$6"
visualize_specific_tgfm_examples_dir="$7"


if false; then
python3 visualize_specific_tgfm_examples_preprocessing.py $specific_examples_input_file $tgfm_input_summary_file $tgfm_results_dir $tgfm_organized_results_dir $gtex_susie_gene_models_dir $gene_annotation_file $visualize_specific_tgfm_examples_dir
fi


Rscript visualize_specific_tgfm_examples.R $specific_examples_input_file $visualize_specific_tgfm_examples_dir