#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)








gtex_pseudotissue_file="$1"
pb_cell_type_file="$2"
preprocessed_tgfm_sldsc_data_dir="$3"
gene_type="$4"
gtex_susie_gene_models_dir="$5"
sc_pbmc_gene_models_dir="$6"
source ~/.bash_profile



tissue_version="no_testis"
python3 organize_gene_ld_scores_for_sc_tgfm_sldsc.py $gtex_pseudotissue_file $pb_cell_type_file $preprocessed_tgfm_sldsc_data_dir $gene_type $gtex_susie_gene_models_dir $sc_pbmc_gene_models_dir $tissue_version