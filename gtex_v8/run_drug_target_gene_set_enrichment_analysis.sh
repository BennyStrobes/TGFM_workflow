#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-7:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=15GB                         # Memory total in MiB (for all cores)


tgfm_results_dir="$1"
gene_type="$2"
trait_name_file="$3"
gene_annotation_file="$4"
gtex_susie_gene_models_dir="$5"
preprocessed_tgfm_data_dir="$6"
drug_target_gene_list_file="$7"
drug_target_gene_set_enrichment_dir="$8"


python3 run_drug_target_gene_set_enrichment_analysis.py $tgfm_results_dir $gene_type $trait_name_file $gene_annotation_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $drug_target_gene_list_file $drug_target_gene_set_enrichment_dir


if false; then

module load R/3.5.1

Rscript visualize_drug_target_gene_set_enrichment_analysis.R $drug_target_gene_set_enrichment_dir

python3 run_mendelian_gene_set_enrichment_analysis.py $tgfm_results_dir $gene_type $trait_name_file $gene_annotation_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $drug_target_gene_list_file $drug_target_gene_set_enrichment_dir


Rscript visualize_mendelian_gene_set_enrichment_analysis.R $drug_target_gene_set_enrichment_dir
fi