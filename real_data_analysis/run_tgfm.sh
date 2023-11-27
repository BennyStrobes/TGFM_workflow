#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)



trait_name="$1"
ukkbb_window_summary_file="$2"
gtex_pseudotissue_file="$3"
preprocessed_tgfm_data_dir="$4"
tgfm_sldsc_results_dir="$5"
samp_size="$6"
gene_type="$7"
tgfm_results_dir="$8"
job_number="$9"
num_jobs="${10}"




source ~/.bash_profile
module load R/4.0.1

echo $trait_name
echo $job_number



python3 run_tgfm.py $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $tgfm_sldsc_results_dir $samp_size $gene_type $tgfm_results_dir $job_number $num_jobs