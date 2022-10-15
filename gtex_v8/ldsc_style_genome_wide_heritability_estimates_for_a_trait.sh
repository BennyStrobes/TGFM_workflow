#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





trait_name="$1"
ukkbb_window_summary_file="$2"
tissue_name_file="$3"
preprocessed_tgfm_data_dir="$4"
output_stem="$5"


module load gcc/6.2.0
module load python/3.6.0
source /n/groups/price/ben/environments/tensor_flow_cpu/bin/activate


python3 ldsc_style_genome_wide_heritability_estimates_for_a_trait.py $trait_name $ukkbb_window_summary_file $tissue_name_file $preprocessed_tgfm_data_dir $output_stem
