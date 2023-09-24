#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4GB                         # Memory total in MiB (for all cores)




pseudotissue_name="$1"
gtex_pseudotissue_gene_model_input_dir="$2"
gtex_susie_gene_models_dir="$3"

source ~/.bash_profile

echo $pseudotissue_name


python3 organize_susie_gene_model_results_in_a_single_pseudotissue.py $pseudotissue_name $gtex_pseudotissue_gene_model_input_dir $gtex_susie_gene_models_dir