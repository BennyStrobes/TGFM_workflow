#!/bin/bash
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)



trait_name="$1"
pseudotissue_gtex_rss_multivariate_twas_dir="$2"
gene_version="$3"
samp_size="$4"





module load gcc/6.2.0
module load python/3.6.0
source /n/groups/price/ben/environments/tensor_flow_cpu/bin/activate




python3 run_ldsc_robust_rss_twas_regression.py $trait_name $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $samp_size