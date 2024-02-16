#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)




tgfm_results_dir="$1"
gene_type="$2"
num_jobs="$3"
trait_names_file="$4"
gtex_pseudotissue_file="$5"
gtex_pseudotissue_category_file="$6"
processed_tgfm_input_stem="$7"
ukbb_preprocessed_for_genome_wide_susie_dir="$8"
tgfm_organized_results_dir="${9}"
gene_annotation_file="${10}"
tgfm_gtex_organized_results_dir="${11}"


source ~/.bash_profile

python3 organize_sc_tgfm_results_across_parallel_runs.py $tgfm_results_dir $gene_type $num_jobs $trait_names_file $gtex_pseudotissue_file $gtex_pseudotissue_category_file $processed_tgfm_input_stem $ukbb_preprocessed_for_genome_wide_susie_dir $tgfm_organized_results_dir $gene_annotation_file $tgfm_gtex_organized_results_dir