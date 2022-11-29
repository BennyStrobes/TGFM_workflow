#!/bin/bash
#SBATCH -c 20                               # Request one core
#SBATCH -t 0-60:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30G                         # Memory total in MiB (for all cores)


module load gcc/6.2.0
module load python/3.6.0
source /n/groups/price/ben/environments/tensor_flow_cpu/bin/activate



trait_name="$1"
ukkbb_window_summary_file="$2"
gtex_pseudotissue_file="$3"
preprocessed_tgfm_data_dir="$4"
standardize_expression_boolean="$5"
output_root="$6"
gene_type="$7"


python3 rss_likelihood_parallel_genome_wide_heritability_estimates_for_a_trait.py $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $standardize_expression_boolean $output_root $gene_type
