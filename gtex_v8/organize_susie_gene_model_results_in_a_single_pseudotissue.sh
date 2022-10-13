#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




pseudotissue_name="$1"
gtex_pseudotissue_gene_model_input_dir="$2"
gtex_susie_gene_models_dir="$3"

if false; then
source ~/.bash_profile
fi


python3 organize_susie_gene_model_results_in_a_single_pseudotissue.py $pseudotissue_name $gtex_pseudotissue_gene_model_input_dir $gtex_susie_gene_models_dir