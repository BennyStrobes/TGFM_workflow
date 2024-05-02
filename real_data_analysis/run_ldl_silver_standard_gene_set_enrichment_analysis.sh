#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=15GB                         # Memory total in MiB (for all cores)





tgfm_organized_results_dir="$1"
gtex_susie_gene_models_dir="$2"
gene_annotation_file="$3"
ldl_silver_standard_gene_set_file="$4"
ldl_silver_standard_gene_set_enrichment_dir="$5"

python3 run_ldl_silver_standard_gene_set_enrichment_analysis.py $tgfm_organized_results_dir $gtex_susie_gene_models_dir $gene_annotation_file $ldl_silver_standard_gene_set_file $ldl_silver_standard_gene_set_enrichment_dir


python3 run_ldl_silver_standard_gene_tissue_set_enrichment_analysis.py $tgfm_organized_results_dir $gtex_susie_gene_models_dir $gene_annotation_file $ldl_silver_standard_gene_set_file $ldl_silver_standard_gene_set_enrichment_dir
