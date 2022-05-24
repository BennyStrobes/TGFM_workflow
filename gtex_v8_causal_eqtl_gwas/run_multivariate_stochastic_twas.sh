#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-80:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)




trait_name="$1"
ukbb_genome_wide_susie_organized_results_dir="$2"
gtex_pseudotissue_file="$3"
pseudotissue_gtex_susie_pmces_fusion_weights_dir="$4"
gtex_fusion_processed_intermediate_data="$5"
samp_size="$6"
pseudotissue_gtex_stochastic_multivariate_twas_dir="$7"
chrom_num="$8"

echo $chrom_num


module load R/4.0.1


Rscript run_multivariate_stochastic_twas.R $chrom_num $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_stochastic_multivariate_twas_dir

