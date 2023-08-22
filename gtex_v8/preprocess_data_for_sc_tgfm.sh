#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-14:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





ukkbb_window_summary_file="$1"
hapmap3_snpid_file="$2"
gtex_pseudotissue_file="$3"
gtex_susie_gene_models_dir="$4"
preprocessed_tgfm_data_dir="$5"
job_number="$6"
num_jobs="$7"
gene_type="$8"
annotation_dir="$9"
pb_cell_type_file="${10}"
sc_pbmc_susie_gene_models_dir="${11}"

source ~/.bash_profile

echo $job_number

date

python3 preprocess_data_for_sc_tgfm.py $ukkbb_window_summary_file $hapmap3_snpid_file $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $job_number $num_jobs $gene_type $annotation_dir $pb_cell_type_file $sc_pbmc_susie_gene_models_dir

date