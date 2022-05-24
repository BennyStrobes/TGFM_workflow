#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)



trait_name="$1"
ukbb_genome_wide_susie_organized_results_dir="$2"
gtex_tissue_file="$3"
gtex_fusion_weights_dir="$4"
kg_genotype_dir="$5"
gtex_fusion_associations_dir="$6"
samp_size="$7"
gtex_fusion_multivariate_associations_dir="$8"


if false; then
for chrom_num in {1..22}; do 
	sbatch run_multivariate_twas_from_summary_stats_with_susie_for_single_trait_and_chromosome.sh $chrom_num $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_tissue_file $gtex_fusion_weights_dir $kg_genotype_dir $gtex_fusion_associations_dir $samp_size $gtex_fusion_multivariate_associations_dir
done
fi

echo $trait_name
module load R/4.0.1
Rscript organize_multivariate_twas_results.R $trait_name $gtex_fusion_multivariate_associations_dir $ukbb_genome_wide_susie_organized_results_dir $gtex_tissue_file
