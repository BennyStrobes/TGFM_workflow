#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-14:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)





ukkbb_window_summary_file="$1"
gtex_pseudotissue_file="$2"
gtex_susie_gene_models_dir="$3"
preprocessed_tgfm_data_dir="$4"
job_number="$5"
num_jobs="$6"
gene_type="$7"
pb_cell_type_file="${8}"
sc_pbmc_susie_gene_models_dir="${9}"

source ~/.bash_profile

echo $job_number

date

python3 preprocess_data_for_sc_tgfm.py $ukkbb_window_summary_file $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $job_number $num_jobs $gene_type $pb_cell_type_file $sc_pbmc_susie_gene_models_dir

date