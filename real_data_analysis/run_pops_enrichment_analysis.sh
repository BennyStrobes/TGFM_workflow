#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=15GB                         # Memory total in MiB (for all cores)


tgfm_results_dir="$1"
gene_type="$2"
trait_names_file="$3"
gtex_susie_gene_models_dir="$4"
preprocessed_tgfm_data_dir="$5"
pops_results_summary_file="$6"
pops_enrichment_dir="$7"
tgfm_organized_results_dir="$8"

source ~/.bash_profile
python3 run_pops_enrichment_analysis.py $tgfm_results_dir $gene_type $trait_names_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $pops_results_summary_file $pops_enrichment_dir $tgfm_organized_results_dir


if false; then
module load R/3.5.1
Rscript visualize_pops_enrichment.R $pops_enrichment_dir
fi