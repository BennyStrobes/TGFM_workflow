#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=15GB                         # Memory total in MiB (for all cores)



tgfm_results_dir="$1"
traits_file="$2"
gtex_susie_gene_models_dir="$3"
preprocessed_tgfm_data_dir="$4"
tgfm_organized_results_dir="$5"
non_disease_specific_gene_sets_file="$6"
em_gene_set_file="$7"
non_disease_specific_gene_set_enrichment_dir="$8"



python3 run_non_disease_specific_gene_set_enrichment_analysis.py $tgfm_results_dir $traits_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $tgfm_organized_results_dir $non_disease_specific_gene_sets_file $em_gene_set_file $non_disease_specific_gene_set_enrichment_dir